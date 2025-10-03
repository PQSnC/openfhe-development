// rowdot_bfv.cpp
// Build: CMake / C++17 이상
// 목적: 512차원 벡터 × (8129명) 내적을 행 기준 배치로 수행.
//       - 각 행(512)마다 8129명의 값을 하나의 암호문으로 패킹
//       - 샘플 벡터의 스칼라를 모든 슬롯에 복제한 평문과 ct-pt 곱: 512회
//       - 곱 결과들을 ct-ct 덧셈으로 누적: 511회
//       - 누적 최종 암호문 복호화 시간 측정
//
// 사용 예:
//   ./rowdot_bfv 8192 512 8129 3 65537 1 60
//         N     vec   cols runs    t  dep dcrtBits
//
// 주의:
//  - num_cols(8129) ≤ N(8192) 이어야 함. (여유 슬롯은 0으로 채워짐)
//  - t=65537은 full batching 조건(t ≡ 1 (mod 2N))을 만족(N=8192 → 2N=16384, 65536=4*16384)

#include "openfhe.h"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <ios>
#include <iomanip>

using namespace lbcrypto;

template <class F>
static double time_ms(F&& fn) {
    auto t0 = std::chrono::high_resolution_clock::now();
    fn();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

static bool isPow2(size_t x){ return x && ((x & (x-1))==0); }

int main(int argc, char** argv) {
    // args: N vec_len num_cols runs [t] [depth] [dcrtBits]
    uint32_t N        = (argc > 1)? std::stoul(argv[1])  : 8192;   // ring dimension
    size_t   vec_len  = (argc > 2)? std::stoull(argv[2]) : 512;    // 차원(행 수)
    size_t   num_cols = (argc > 3)? std::stoull(argv[3]) : 8129;   // 사람 수(열 수 = 슬롯 사용 수)
    int      runs     = (argc > 4)? std::stoi(argv[4])   : 1;       // 반복 횟수(평균 내기)
    uint64_t t        = (argc > 5)? std::stoull(argv[5]) : 65537;  // 평문 모듈러스
    int      depth    = (argc > 6)? std::stoi(argv[6])   : 1;       // 필요 곱 깊이(여기선 1)
    int      dcrtBits = (argc > 7)? std::stoi(argv[7])   : 60;      // CRT 소수 비트

    if (!isPow2(N)) {
        std::cerr << "[error] N must be a power of two.\n";
        return 1;
    }
    if (num_cols > N) {
        std::cerr << "[error] num_cols must be ≤ N (slots). Given num_cols="<<num_cols<<", N="<<N<<"\n";
        return 1;
    }
    if (vec_len < 1) {
        std::cerr << "[error] vec_len must be ≥ 1.\n";
        return 1;
    }

    // ---- BFV 컨텍스트 생성 (수동 파라미터) ----
    CCParams<CryptoContextBFVRNS> ps;
    ps.SetSecurityLevel(HEStd_NotSet);    // 수동 설정
    ps.SetRingDim(N);
    ps.SetBatchSize(0);                   // full batching (num_cols가 2의 거듭제곱이 아니므로 0 권장)
    ps.SetPlaintextModulus(t);
    ps.SetMultiplicativeDepth(depth);     // ct*pt 1회면 충분
    ps.SetScalingModSize(dcrtBits);

    auto cc = GenCryptoContext(ps);
    cc->Enable(PKE);
    cc->Enable(LEVELEDSHE);

    // 파라미터 요약
    auto N_eff = cc->GetRingDimension();
    auto ep = cc->GetCryptoParameters()->GetElementParams();
    using ParamsT = DCRTPoly::Params;
    auto dcrtParams = std::dynamic_pointer_cast<ParamsT>(ep);
    size_t numPrimes = dcrtParams ? dcrtParams->GetParams().size() : 0;
    double logQ = 0.0;
    if (dcrtParams){
        for (const auto& p : dcrtParams->GetParams())
            logQ += static_cast<double>(p->GetModulus().GetMSB());
    }

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "=== BFV Row-Batched Dot Timing ===\n"
              << "N="<<N<<" (eff "<<N_eff<<"), dims(vec_len)="<<vec_len
              << ", persons(num_cols)="<<num_cols<<", runs="<<runs<<"\n"
              << "t="<<t<<", depth="<<depth<<", dcrtBits="<<dcrtBits
              << " -> numCRTprimes="<<numPrimes<<", approx log2(Q)="<<std::setprecision(1)<<logQ<<" bits\n\n";

    // ---- 키 생성 ----
    auto keys = cc->KeyGen();
    if (!keys.good()) {
        std::cerr << "KeyGen failed.\n";
        return 1;
    }

    // ---- DB 행 기준 배치: 각 행 i에 대해 (num_cols명) 값을 하나의 암호문에 패킹 ----
    // 데이터는 데모용으로 [-1000,1000] 난수 사용 (실데이터로 교체 가능)
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int64_t> dist(-1000, 1000);

    std::vector<Ciphertext<DCRTPoly>> ctRows(vec_len);
    {
        // 행 i: 길이 num_cols의 벡터를 패킹 → 평문화 → 암호화
        for (size_t i=0;i<vec_len;++i) {
            std::vector<int64_t> row(num_cols);
            for (size_t j=0;j<num_cols;++j) row[j] = dist(rng); // DB 값 (32bit 정수라면 범위를 조절)
            Plaintext ptRow = cc->MakePackedPlaintext(row);
            ctRows[i] = cc->Encrypt(keys.publicKey, ptRow);
        }
    }

    // ---- 워밍업 (JIT / 메모리 캐시 영향 제거) ----
    {
        // 샘플 스칼라 벡터 생성
        std::vector<int64_t> q(vec_len);
        for (size_t i=0;i<vec_len;++i) q[i] = dist(rng);

        // 첫 행만 ct*pt, 그리고 덧셈/복호화를 한 번 수행
        Plaintext ptScalar0 = cc->MakePackedPlaintext(std::vector<int64_t>(num_cols, q[0]));
        auto warmMul = cc->EvalMult(ctRows[0], ptScalar0);
        auto warmAcc = warmMul;
        if (vec_len > 1) warmAcc = cc->EvalAdd(warmAcc, warmMul);
        Plaintext warmOut;
        cc->Decrypt(keys.secretKey, warmAcc, &warmOut);
    }

    // ---- 본 측정 ----
    double sum_ms_mul = 0.0;   // 모든 run에서 512번 ct*pt의 총 시간
    double sum_ms_add = 0.0;   // 모든 run에서 511번 ct+ct의 총 시간
    double sum_ms_dec = 0.0;   // 모든 run에서 최종 복호(1회)의 총 시간
    double sum_ms_total = 0.0; // end-to-end(곱+덧셈+복호)의 총 시간

    for (int r=0; r<runs; ++r) {
        // 샘플 512차원 벡터 (평문). 각 스칼라 q[i]는 모든 슬롯에 복제되어 ct*pt에 사용됨.
        std::vector<int64_t> q(vec_len);
        for (size_t i=0;i<vec_len;++i) q[i] = dist(rng);

        // ct*pt 512회
        std::vector<Ciphertext<DCRTPoly>> partial(vec_len);
        double ms_mul = time_ms([&]{
            for (size_t i=0;i<vec_len;++i) {
                Plaintext ptScalar = cc->MakePackedPlaintext(std::vector<int64_t>(num_cols, q[i]));
                partial[i] = cc->EvalMult(ctRows[i], ptScalar);
            }
        });

        // ct+ct 511회 (부분합 누적 → 최종 암호문 acc)
        Ciphertext<DCRTPoly> acc = partial[0];
        double ms_add = time_ms([&]{
            for (size_t i=1;i<vec_len;++i) {
                acc = cc->EvalAdd(acc, partial[i]);
            }
        });

        // 최종 복호 1회
        double ms_dec = 0.0;
        Plaintext out;
        ms_dec += time_ms([&]{ cc->Decrypt(keys.secretKey, acc, &out); });

        sum_ms_mul   += ms_mul;
        sum_ms_add   += ms_add;
        sum_ms_dec   += ms_dec;
        sum_ms_total += (ms_mul + ms_add + ms_dec);

        // (옵션) 결과 sanity check: 앞 8명만 출력 (모듈러 t 값)
        if (r == runs-1) {
            out->SetLength(num_cols);
            auto vec = out->GetPackedValue();
            std::cout << "[Sanity] first 8 persons (mod t): ";
            for (int j=0; j<8 && j<(int)num_cols; ++j) std::cout << vec[j] << (j==7?'\n':' ');
        }
    }

    // ---- 결과 출력 ----
    const double runs_d = static_cast<double>(runs);
    const double per_mul = sum_ms_mul / (runs_d * std::max<size_t>(1, vec_len));         // 멀티플라이 1회당 평균
    const double per_add = sum_ms_add / (runs_d * std::max<size_t>(1, vec_len>0?vec_len-1:0)); // 덧셈 1회당 평균

    std::cout << std::setprecision(3);
    std::cout << "\n[Workload per run]\n"
              << "  ct*pt multiplies : " << vec_len  << " (row-wise)\n"
              << "  ct+ct adds       : " << (vec_len>0?vec_len-1:0) << "\n"
              << "  decrypts         : 1\n\n";

    std::cout << "[Timing]\n"
              << "  PlainMult  total per run : " << (sum_ms_mul / runs_d) << " ms"
              << "  (avg per ct*pt: " << per_mul << " ms)\n"
              << "  CtAdd       total per run : " << (sum_ms_add / runs_d) << " ms"
              << "  (avg per add : " << per_add << " ms)\n"
              << "  Decrypt     per run       : " << (sum_ms_dec / runs_d) << " ms\n"
              << "  TOTAL       per run       : " << (sum_ms_total / runs_d) << " ms\n";

    std::cout << std::endl;
    return 0;
}
