// test_qparam_patched.cpp
#include "openfhe.h"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <ios>

using namespace lbcrypto;

template <class F>
inline double ms(F&& f){
    auto t0 = std::chrono::high_resolution_clock::now();
    f();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count()/1000.0;
}

template <class F>
static double avg_ms(F&& fn, int runs) {
    using clock = std::chrono::high_resolution_clock;
    double total_ms = 0.0;
    for (int i = 0; i < runs; ++i) {
        auto t0 = clock::now();
        fn();
        auto t1 = clock::now();
        total_ms += std::chrono::duration<double, std::milli>(t1 - t0).count();
    }
    return total_ms / runs;
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
    // args: N vec_len num_cols runs [t] [depth] [dcrtBits] [ks] [mul] [measureRot] [autoMode] [constK]
    uint32_t N          = (argc > 1)? std::stoul(argv[1])  : 2048;
    size_t   vec_len    = (argc > 2)? std::stoull(argv[2]) : 512;
    size_t   num_cols   = (argc > 3)? std::stoull(argv[3]) : 2048;
    int      runs       = (argc > 4)? std::stoi(argv[4])   : 1;
    uint64_t t          = (argc > 5)? std::stoull(argv[5]) : 65537;
    int      depth      = (argc > 6)? std::stoi(argv[6])   : 1;
    int      dcrtBits   = (argc > 7)? std::stoi(argv[7])   : 55;
    int      ksArg      = (argc > 8)? std::stoi(argv[8])   : 0;
    int      mulArg     = (argc > 9)? std::stoi(argv[9])   : 0;
    int      measureRot = (argc >10)? std::stoi(argv[10])  : 0;
    int      autoMode   = (argc >11)? std::stoi(argv[11])  : 0;
    long long kConst    = (argc >12)? std::stoll(argv[12]) : 5;

    CryptoContext<DCRTPoly> cc;
    size_t requested_vec_len = vec_len;

    if (autoMode) {
        // === 첫 번째 코드와 같은 자동 파라미터 ===
        CCParams<CryptoContextBFVRNS> params;
        params.SetSecurityLevel(HEStd_128_classic);
        params.SetPlaintextModulus(t);
        params.SetMultiplicativeDepth(depth);
        cc = GenCryptoContext(params);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);

        // 실제 선택된 N/슬롯
        N = cc->GetRingDimension();
        size_t nslots = N;
        // 첫 코드와 동일하게 슬롯 수를 제한(최대 4096)
        vec_len = std::min(nslots, (isPow2(requested_vec_len) && requested_vec_len>0) ? requested_vec_len : size_t(4096));
    } else {
        // === 수동 파라미터 ===
        if (!isPow2(N)) { std::cerr << "[error] N must be a power of two.\n"; return 1; }
        if (vec_len > N) { std::cerr << "[warn] vec_len > N — clamping to N\n"; vec_len = N; }
        size_t batchSize = isPow2(vec_len) ? vec_len : 0;
        if (batchSize==0 && !isPow2(vec_len))
            std::cerr << "[info] vec_len="<<vec_len<<" not power-of-two; using BatchSize=full(0). Extra slots will be zero.\n";
        if (((t-1) % (2ULL*N)) != 0ULL)
            std::cerr << "[warn] t may not satisfy t ≡ 1 (mod 2N). Full batching may be disabled.\n";

        CCParams<CryptoContextBFVRNS> ps;
        ps.SetSecurityLevel(HEStd_NotSet);
        ps.SetRingDim(N);
        ps.SetBatchSize(batchSize);
        ps.SetPlaintextModulus(t);
        ps.SetMultiplicativeDepth(depth);
        ps.SetScalingModSize(dcrtBits);
        ps.SetKeySwitchTechnique(parseKS(ksArg));
        ps.SetMultiplicationTechnique(parseMT(mulArg));

        cc = GenCryptoContext(ps);
        cc->Enable(PKE);
        cc->Enable(KEYSWITCH);
        cc->Enable(LEVELEDSHE);
        cc->Enable(ADVANCEDSHE);
    }

    // ---- 실제 채택된 N, Q 요약 출력 ----
    auto N_eff = cc->GetRingDimension();
    auto ep = cc->GetCryptoParameters()->GetElementParams();
    using ParamsT = DCRTPoly::Params;
    auto dcrtParams = std::dynamic_pointer_cast<ParamsT>(ep);

    size_t numPrimes = dcrtParams ? dcrtParams->GetParams().size() : 0;
    double logQbits = 0.0;
    if (dcrtParams){
        for (const auto& p : dcrtParams->GetParams())
            logQbits += static_cast<double>(p->GetModulus().GetMSB()); // ≈ log2(modulus)
    }

    std::cout.setf(std::ios::fixed);
    std::cout << "=== OpenFHE BFV Inner-Product Timing ==="
              << (autoMode ? " [AUTO]" : " [MANUAL]") << "\n";
    std::cout << "requested N=" << (autoMode ? 0 : N)
              << ", effective N=" << N_eff
              << ", slots(vec_len)=" << vec_len
              << ", persons=" << num_cols
              << ", runs=" << runs << "\n";
    std::cout << "t=" << t << ", depth=" << depth << ", dcrtBits=" << dcrtBits
              << ", ks=" << ((ksArg==1)?"HYBRID":"BV")
              << ", mul=" << ((mulArg==0)?"HPS":(mulArg==1)?"HPSPOVERQ":(mulArg==2)?"HPSPOVERQLEVELED":"BEHZ")
              << " -> numCRTprimes=" << numPrimes
              << ", approx log2(Q)=" << std::setprecision(1) << logQbits << " bits\n\n";

    // ---- 키 & 회전키 ----
    auto keys = cc->KeyGen();
    if (!keys.good()) { std::cerr << "KeyGen failed\n"; return 1; }

    std::vector<int32_t> rotIdx;
    for (size_t step=1; step<vec_len; step<<=1) rotIdx.push_back((int32_t)step);
    if (!rotIdx.empty())
        cc->EvalAtIndexKeyGen(keys.secretKey, rotIdx);

    // ---- (A) 첫 번째 코드와 동일한 마이크로 벤치 3종 ----
    {
        std::vector<int64_t> a(vec_len), b(vec_len);
        for (size_t i = 0; i < vec_len; ++i) {
            a[i] = static_cast<int64_t>((i % 10) + 1);
            b[i] = static_cast<int64_t>((i % 7) + 1);
        }
        Plaintext pta = cc->MakePackedPlaintext(a);
        Plaintext ptb = cc->MakePackedPlaintext(b);
        auto ctA = cc->Encrypt(keys.publicKey, pta);
        auto ctB = cc->Encrypt(keys.publicKey, ptb);
        Plaintext ptConst = cc->MakePackedPlaintext(std::vector<int64_t>(vec_len, (int64_t)kConst));

        // 워밍업
        (void)cc->EvalAdd(ctA, ctB);
        (void)cc->EvalMult(ctA, ptConst);
        Plaintext tmp; cc->Decrypt(keys.secretKey, ctA, &tmp);

        double add_ms_once = avg_ms([&]{ auto ct = cc->EvalAdd(ctA, ctB); (void)ct; }, runs);
        double mul_const_ms = avg_ms([&]{ auto ct = cc->EvalMult(ctA, ptConst); (void)ct; }, runs);
        double dec_ms_once = avg_ms([&]{ Plaintext out; cc->Decrypt(keys.secretKey, ctA, &out); }, runs);

        std::cout << "[Match-First] EvalAdd (ct+ct)            : " << std::setprecision(3) << add_ms_once << " ms\n";
        std::cout << "[Match-First] EvalMult (ct*const via pt) : " << mul_const_ms << " ms\n";
        std::cout << "[Match-First] Decrypt                    : " << dec_ms_once << " ms\n\n";
    }

    // ---- (B) 원래 하시던 inner-product 스타일 측정 ----
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int64_t> dist(-1000, 1000);

    std::vector<int64_t> q(vec_len); for (size_t i=0;i<vec_len;++i) q[i]=dist(rng);
    Plaintext ptQuery = cc->MakePackedPlaintext(q);

    double sum_mul=0.0, sum_rot=0.0, sum_add=0.0, sum_dec=0.0;

    // 워밍업
    {
        std::vector<int64_t> v(vec_len);
        for(size_t i=0;i<vec_len;++i) v[i]=dist(rng);
        auto ct = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(v));
        auto had = cc->EvalMult(ct, ptQuery);
        if (!rotIdx.empty()) {
            auto r1  = cc->EvalAtIndex(had, 1);
            (void)cc->EvalAdd(had, r1);
        }
        Plaintext d; cc->Decrypt(keys.secretKey, had, &d);
    }

    for(size_t col=0; col<num_cols; ++col){
        std::vector<int64_t> person(vec_len);
        for(size_t i=0;i<vec_len;++i) person[i]=dist(rng);
        auto ctPerson = cc->Encrypt(keys.publicKey, cc->MakePackedPlaintext(person));

        for(int r=0; r<runs; ++r){
            // (1) ct * pt (성분곱)
            Ciphertext<DCRTPoly> ctHad;
            sum_mul += ms([&]{ ctHad = cc->EvalMult(ctPerson, ptQuery); });

            // (2) 슬롯 합산: 회전 + 덧셈
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
    double add_per_call_ms = sum_add / denom / std::log2((double)vec_len);

    std::cout << std::setprecision(3);
    std::cout << "[PlainMult] ct*pt (slotwise)  avg: " << (sum_mul/denom) << " ms\n";
    if (measureRot) {
        std::cout << "[Rotate]    EvalAtIndex       avg: " << (sum_rot/denom) << " ms\n";
        std::cout << "[Rot+Add]   slot-sum total    avg: " << ((sum_rot+sum_add)/denom) << " ms\n";
    }
    std::cout << "[CtAdd]     slot-sum (adds)   avg: " << (sum_add/denom) << " ms  "
              << "(≈ per add " << add_per_call_ms << " ms)\n";
    std::cout << "[Decrypt]   final scalar      avg: " << (sum_dec/denom) << " ms\n";
    if (measureRot) {
        std::cout << "[Total]     end-to-end        avg: "
                  << ((sum_mul+sum_rot+sum_add+sum_dec)/denom) << " ms\n";
    }
    return 0;
}
