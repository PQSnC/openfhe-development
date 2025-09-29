// src/pke/examples/bfv-bench.cpp
#include "openfhe.h"
#include <chrono>
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

int main(int argc, char** argv) {
    // 인자: ringDim vec_len num_cols runs [k]
    //  - ringDim: 링 차원 N (2의 거듭제곱 권장; 예 2048)
    //  - vec_len:  슬롯(배치) 수 (예 512)
    //  - num_cols: 사람 수(열 개수, 암호문 개수; 예 2048)
    //  - runs:     같은 열에 대해 반복 측정 횟수
    //  - k:        (이전 호환용) 사용하지 않음
    uint32_t ringDim = (argc > 1)? std::stoul(argv[1]) : 2048;
    size_t   vec_len = (argc > 2)? std::stoull(argv[2]) : 512;
    size_t   num_cols= (argc > 3)? std::stoull(argv[3]) : 2048;
    int      runs    = (argc > 4)? std::stoi(argv[4])   : 1;
    // int64_t  k       = (argc > 5)? std::stoll(argv[5])  : 7; // <- 내적 시나리오에서는 사용 안 함

    // === BFV 컨텍스트 설정 ===
    CCParams<CryptoContextBFVRNS> ps;
    ps.SetRingDim(ringDim);
    ps.SetBatchSize(vec_len);             // 슬롯 수 = 512 등
    ps.SetPlaintextModulus(65537);        // 배칭용: 65537 ≡ 1 (mod 2N)
    ps.SetMultiplicativeDepth(1);         // ct×pt, 덧셈만 필요
    ps.SetSecurityLevel(SecurityLevel::HEStd_NotSet); // 벤치용(수동 파라미터)

    auto cc = GenCryptoContext(ps);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);   // 회전(AtIndex)에 필요
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE); // 일부 환경에서 회전/합산에 필요

    // 키 생성
    auto keys = cc->KeyGen();
    if (!keys.good()) {
        std::cerr << "KeyGen failed\n";
        return 1;
    }

    // 회전 키(1,2,4,...) 생성: 슬롯 합산(rotate-sum)에 필요
    std::vector<int32_t> rotIdx;
    for (size_t step = 1; step < vec_len; step <<= 1)
        rotIdx.push_back(static_cast<int32_t>(step));
    cc->EvalAtIndexKeyGen(keys.secretKey, rotIdx);

    // ===== 입력 쿼리(평문 1×vec_len) 준비 =====
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int64_t> dist(-1000, 1000);

    std::vector<int64_t> q(vec_len);
    for (size_t i = 0; i < vec_len; ++i) q[i] = dist(rng);
    Plaintext ptQuery = cc->MakePackedPlaintext(q);

    // ===== 누적 시간 =====
    double sum_mul = 0.0;     // ct × pt (슬롯별 곱)
    double sum_add = 0.0;     // rotate-sum에서의 "덧셈" 시간(회전 제외)
    double sum_dec = 0.0;     // 최종(합산 결과) 복호화

    // ===== 워밍업 =====
    {
        std::vector<int64_t> v(vec_len);
        for (size_t i = 0; i < vec_len; ++i) v[i] = dist(rng);
        auto pt = cc->MakePackedPlaintext(v);
        auto ct = cc->Encrypt(keys.publicKey, pt);
        auto tmpMul = cc->EvalMult(ct, ptQuery); // ct×pt
        auto rot = cc->EvalAtIndex(tmpMul, 1);
        (void)cc->EvalAdd(tmpMul, rot);
        Plaintext d; cc->Decrypt(keys.secretKey, tmpMul, &d);
    }

    // ===== 본 측정: 사람(열)마다 내적 계산 =====
    for (size_t col = 0; col < num_cols; ++col) {
        // 이 사람의 512차원 벡터(평문) 생성 & 암호화 -> 하나의 배치 암호문
        std::vector<int64_t> person(vec_len);
        for (size_t i = 0; i < vec_len; ++i) person[i] = dist(rng);
        Plaintext pPerson = cc->MakePackedPlaintext(person);
        auto ctPerson = cc->Encrypt(keys.publicKey, pPerson);

        for (int r = 0; r < runs; ++r) {
            // (1) ct × pt (슬롯별 곱) : 내적의 Hadamard 곱 단계
            Ciphertext<DCRTPoly> ctHad;
            double t_mul = ms([&]{ ctHad = cc->EvalMult(ctPerson, ptQuery); });
            sum_mul += t_mul;

            // (2) 슬롯 합산: 회전(+회전키 필요) 후 덧셈 누적
            //     요구사항상 "덧셈 시간만" 계측하고, 회전 시간은 제외
            Ciphertext<DCRTPoly> acc = ctHad;
            double t_add_total = 0.0;
            for (size_t step = 1; step < vec_len; step <<= 1) {
                auto rotated = cc->EvalAtIndex(acc, static_cast<int32_t>(step)); // 회전(비계측)
                t_add_total += ms([&]{ acc = cc->EvalAdd(acc, rotated); });      // 덧셈만 계측
            }
            sum_add += t_add_total;

            // (3) 복호화 시간: acc는 모든 슬롯이 내적 결과(스칼라)를 담음
            Plaintext d;
            double t_dec = ms([&]{ cc->Decrypt(keys.secretKey, acc, &d); });
            sum_dec += t_dec;

            // 필요 시 정확성 확인(첫 열에서만 샘플)
            // if (col == 0 && r == 0) {
            //     d->SetLength(1);
            //     auto vals = d->GetPackedValue();
            //     std::cout << "[sample] dot=" << vals[0] << "\n";
            // }
        }
    }

    double denom = static_cast<double>(num_cols * runs);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "=== OpenFHE BFV Inner-Product Timing ===\n";
    std::cout << "ringDim(n)     : " << ringDim << "\n";
    std::cout << "slots(vec_len) : " << vec_len << "\n";
    std::cout << "num persons    : " << num_cols << "  (암호문 개수)\n";
    std::cout << "runs per person: " << runs << "\n\n";
    std::cout << "[PlainMult] ct * pt (slotwise)  avg: " << (sum_mul/denom) << " ms\n";
    std::cout << "[CtAdd]     slot-sum (adds only) avg: " << (sum_add/denom) << " ms\n";
    std::cout << "[Decrypt]   final scalar         avg: " << (sum_dec/denom) << " ms\n";
    return 0;
}
