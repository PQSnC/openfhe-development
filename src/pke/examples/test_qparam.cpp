#include "openfhe.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace lbcrypto;

template <class F>
inline double ms(F&& f){
    auto t0 = std::chrono::high_resolution_clock::now();
    f();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1000.0;
}

static KeySwitchTechnique parseKS(int v){ return (v==1)? HYBRID : BV; }
static MultiplicationTechnique parseMT(int v){
    switch(v){
        case 1: return HPSPOVERQ;
        case 2: return HPSPOVERQLEVELED;
        case 3: return BEHZ;
        default: return HPS;
    }
}

static bool isPow2(size_t x){ return x && ((x & (x-1))==0); }

int main(int argc, char** argv) {
    // args: N vec_len num_cols runs [t] [depth] [dcrtBits] [ks] [mul] [measureRot]
    uint32_t N        = (argc > 1)? std::stoul(argv[1]) : 2048;
    size_t   vec_len  = (argc > 2)? std::stoull(argv[2]) : 512;
    size_t   num_cols = (argc > 3)? std::stoull(argv[3]) : 2048;
    int      runs     = (argc > 4)? std::stoi(argv[4])   : 1;
    uint64_t t        = (argc > 5)? std::stoull(argv[5]) : 65537;
    int      depth    = (argc > 6)? std::stoi(argv[6])   : 1;
    int      dcrtBits = (argc > 7)? std::stoi(argv[7])   : 55;
    int      ksArg    = (argc > 8)? std::stoi(argv[8])   : 0;
    int      mulArg   = (argc > 9)? std::stoi(argv[9])   : 0;
    int      measureRot = (argc >10)? std::stoi(argv[10]): 0;

    if (!isPow2(N)) { std::cerr << "[error] N must be a power of two.\n"; return 1; }
    if (vec_len > N) { std::cerr << "[warn] vec_len > N — clamping to N\n"; vec_len = N; }
    // BFV: BatchSize는 0(full) 또는 2의 거듭제곱만 허용
    size_t batchSize = isPow2(vec_len) ? vec_len : 0;
    if (batchSize==0 && !isPow2(vec_len))
        std::cerr << "[info] vec_len="<<vec_len<<" not power-of-two; using BatchSize=full(0). Extra slots will be zero.\n";
    // 배치 조건 권고: t ≡ 1 (mod 2N)
    if (((t-1) % (2ULL*N)) != 0ULL)
        std::cerr << "[warn] t may not satisfy t ≡ 1 (mod 2N). Full batching may be disabled.\n";

    // ---- BFV context (수동 파라미터) ----
    CCParams<CryptoContextBFVRNS> ps;
    ps.SetSecurityLevel(HEStd_NotSet);   // 수동
    ps.SetRingDim(N);
    ps.SetBatchSize(batchSize);
    ps.SetPlaintextModulus(t);
    ps.SetMultiplicativeDepth(depth);    // 깊이 -> CRT 소수 개수 유도
    ps.SetScalingModSize(dcrtBits);      // 각 CRT 소수 비트수
    ps.SetKeySwitchTechnique(parseKS(ksArg));
    ps.SetMultiplicationTechnique(parseMT(mulArg));

    auto cc = GenCryptoContext(ps);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);               // rotations
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);

    // ---- 실제 채택된 N, Q 요약 출력 ----
    auto N_eff = cc->GetRingDimension();
    auto ep = cc->GetCryptoParameters()->GetElementParams();
    auto dcrtParams = std::dynamic_pointer_cast<ILDCRTParams<NativeInteger>>(ep);
    size_t numPrimes = dcrtParams ? dcrtParams->GetParams().size() : 0;
    double logQ = 0.0;
    if (dcrtParams){
        for (const auto& p : dcrtParams->GetParams())
            logQ += std::log2(p->GetModulus().ConvertToDouble());
    }
    std::cout << "=== OpenFHE BFV Inner-Product Timing ===\n";
    std::cout << "requested N=" << N << ", effective N=" << N_eff
              << ", slots(vec_len)=" << vec_len << ", persons=" << num_cols
              << ", runs=" << runs << "\n";
    std::cout << "t=" << t << ", depth=" << depth << ", dcrtBits=" << dcrtBits
              << ", ks=" << ((ksArg==1)?"HYBRID":"BV")
              << ", mul=" << ((mulArg==0)?"HPS":(mulArg==1)?"HPSPOVERQ":(mulArg==2)?"HPSPOVERQLEVELED":"BEHZ")
              << " -> numCRTprimes=" << numPrimes
              << ", approx log2(Q)=" << std::fixed << std::setprecision(1) << logQ << " bits\n\n";

    // ---- 키 & 회전키 ----
    auto keys = cc->KeyGen();
    if (!keys.good()) { std::cerr << "KeyGen failed\n"; return 1; }
    std::vector<int32_t> rotIdx;
    for (size_t step=1; step<vec_len; step<<=1) rotIdx.push_back((int32_t)step);
    cc->EvalAtIndexKeyGen(keys.secretKey, rotIdx);

    // ---- 쿼리 평문 ----
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int64_t> dist(-1000, 1000);
    std::vector<int64_t> q(vec_len); for (size_t i=0;i<vec_len;++i) q[i]=dist(rng);
    Plaintext ptQuery = cc->MakePackedPlaintext(q);

    // ---- 누적 시간 ----
    double sum_mul=0.0, sum_rot=0.0, sum_add=0.0, sum_dec=0.0;

    // 워밍업
    {
        std::vector<int64_t> v(vec_len);
        for(size_t i=0;i<vec_len;++i) v[i]=dist(rng);
        auto ct = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(v));
        auto had = cc->EvalMult(ct, ptQuery);
        auto r1  = cc->EvalAtIndex(had, 1);
        (void)cc->EvalAdd(had, r1);
        Plaintext d; cc->Decrypt(keys.secretKey, had, &d);
    }

    // 본 측정
    for(size_t col=0; col<num_cols; ++col){
        std::vector<int64_t> person(vec_len);
        for(size_t i=0;i<vec_len;++i) person[i]=dist(rng);
        auto ctPerson = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(person));

        for(int r=0; r<runs; ++r){
            // (1) ct * pt (성분곱)
            Ciphertext<DCRTPoly> ctHad;
            sum_mul += ms([&]{ ctHad = cc->EvalMult(ctPerson, ptQuery); });

            // (2) 슬롯 합산: 회전 + 덧셈 (옵션: 회전도 계측)
            Ciphertext<DCRTPoly> acc = ctHad;
            double t_rot_total=0.0, t_add_total=0.0;
            for(size_t step=1; step<vec_len; step<<=1){
                Ciphertext<DCRTPoly> rotated;
                if (measureRot) t_rot_total += ms([&]{ rotated = cc->EvalAtIndex(acc, (int32_t)step); });
                else            rotated = cc->EvalAtIndex(acc, (int32_t)step);
                t_add_total += ms([&]{ acc = cc->EvalAdd(acc, rotated); });
            }
            sum_rot += t_rot_total;
            sum_add += t_add_total;

            // (3) 복호화
            Plaintext d;
            sum_dec += ms([&]{ cc->Decrypt(keys.secretKey, acc, &d); });
        }
    }

    double denom = double(num_cols)*double(runs);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "[PlainMult] ct*pt (slotwise)  avg: " << (sum_mul/denom) << " ms\n";
    if (measureRot) {
        std::cout << "[Rotate]    EvalAtIndex       avg: " << (sum_rot/denom) << " ms\n";
        std::cout << "[Rot+Add]   slot-sum total    avg: " << ((sum_rot+sum_add)/denom) << " ms\n";
    }
    std::cout << "[CtAdd]     slot-sum (adds)   avg: " << (sum_add/denom) << " ms\n";
    std::cout << "[Decrypt]   final scalar      avg: " << (sum_dec/denom) << " ms\n";
    if (measureRot) {
        std::cout << "[Total]     end-to-end        avg: "
                  << ((sum_mul+sum_rot+sum_add+sum_dec)/denom) << " ms\n";
    }
    return 0;
}
