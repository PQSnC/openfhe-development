// bfv_bench.cpp
// Build: CMake(아래 참고) / C++17 이상
#include "openfhe.h"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace lbcrypto;

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

int main() {
    // ------------------------------
    // 1) CryptoContext 자동 파라미터 설정 (BFVrns)
    // ------------------------------
    CCParams<CryptoContextBFVRNS> params;
    params.SetSecurityLevel(HEStd_128_classic);  // 기본 보안강도
    params.SetPlaintextModulus(65537);           // BFV용 대표적인 소수 t
    params.SetMultiplicativeDepth(1);            // ct-pt 곱(상수곱)과 덧셈만 사용할 예정

    CryptoContext<DCRTPoly> cc = GenCryptoContext(params); // 자동 파라미터 생성
    // 필요한 기능 활성화
    cc->Enable(PKE);         // 공개키 암호화/복호화
    cc->Enable(LEVELEDSHE);  // EvalAdd/EvalMult 등 SHE 연산
    // (상수곱은 ct-pt 곱으로 키스위치/재선형키 불필요)

    auto kp = cc->KeyGen();
    if (!kp.good()) {
        std::cerr << "KeyGen failed\n";
        return 1;
    }

    // 환경 파라미터 출력
    auto cryptoParams = cc->GetCryptoParameters();
    const auto cyclo = cryptoParams->GetElementParams()->GetCyclotomicOrder();
    const std::size_t nslots = static_cast<std::size_t>(cyclo / 2); // BFV packable slots
    std::cout << "BFV (auto params)\n"
              << "  t (plain modulus): " << cryptoParams->GetPlaintextModulus() << "\n"
              << "  ring dimension n : " << nslots << " slots\n";

    // ------------------------------
    // 2) 테스트 데이터 준비 (패킹)
    // ------------------------------
    const std::size_t used_slots = std::min<std::size_t>(nslots, 4096); // 너무 크게 잡지 않기
    std::vector<int64_t> a(used_slots), b(used_slots);
    for (std::size_t i = 0; i < used_slots; ++i) {
        a[i] = static_cast<int64_t>((i % 10) + 1);
        b[i] = static_cast<int64_t>((i % 7) + 1);
    }

    Plaintext pta = cc->MakePackedPlaintext(a);
    Plaintext ptb = cc->MakePackedPlaintext(b);

    auto ctA = cc->Encrypt(kp.publicKey, pta);
    auto ctB = cc->Encrypt(kp.publicKey, ptb);

    // 상수곱: BFV에서는 "암호문 × 평문" 곱을 사용 (상수를 모든 슬롯에 채운 평문으로 인코딩)
    // OpenFHE API는 EvalMult(ciphertext, plaintext)를 제공합니다.
    // (ct-pt 곱 오버로드는 문서화된 표준 경로입니다.)
    // 참조: LeveledSHEBase::EvalMult(ct, pt) :contentReference[oaicite:1]{index=1}
    const int64_t kConst = 5;
    Plaintext ptConst = cc->MakePackedPlaintext(std::vector<int64_t>(used_slots, kConst));

    // ------------------------------
    // 3) 워밍업 (메모리/캐시 등 초기 오버헤드 제거)
    // ------------------------------
    auto ctAddWarm = cc->EvalAdd(ctA, ctB);
    auto ctMulWarm = cc->EvalMult(ctA, ptConst);
    Plaintext tmp;
    cc->Decrypt(kp.secretKey, ctA, &tmp);

    // ------------------------------
    // 4) 시간 측정
    // ------------------------------
    const int runs = 100;

    // 암호문 + 암호문
    const double add_ms = avg_ms([&]() {
        auto ct = cc->EvalAdd(ctA, ctB);
        (void)ct;
    }, runs);

    // 암호문 × 상수 (암호문 × 평문)
    const double mul_const_ms = avg_ms([&]() {
        auto ct = cc->EvalMult(ctA, ptConst);
        (void)ct;
    }, runs);

    // 복호화
    const double dec_ms = avg_ms([&]() {
        Plaintext out;
        cc->Decrypt(kp.secretKey, ctA, &out);
    }, runs);

    // ------------------------------
    // 5) 결과 출력
    // ------------------------------
    std::cout << "\nSlots used       : " << used_slots << "\n"
              << "Runs per op      : " << runs << "\n\n"
              << "[Average times]\n"
              << "  EvalAdd (ct+ct)            : " << add_ms << " ms\n"
              << "  EvalMult (ct*const via pt) : " << mul_const_ms << " ms\n"
              << "  Decrypt                    : " << dec_ms << " ms\n";

    // 정합성 체크: 예시로 앞 8개 슬롯만 확인
    Plaintext addDec, mulDec;
    cc->Decrypt(kp.secretKey, cc->EvalAdd(ctA, ctB), &addDec);
    cc->Decrypt(kp.secretKey, cc->EvalMult(ctA, ptConst), &mulDec);
    addDec->SetLength(used_slots);
    mulDec->SetLength(used_slots);

    auto addVec = addDec->GetPackedValue();
    auto mulVec = mulDec->GetPackedValue();
    std::cout << "\n[Sample outputs, first 8 slots]\n  add: ";
    for (int i = 0; i < 8 && i < static_cast<int>(used_slots); ++i) std::cout << addVec[i] << (i==7?'\n':' ');
    std::cout << "  mul: ";
    for (int i = 0; i < 8 && i < static_cast<int>(used_slots); ++i) std::cout << mulVec[i] << (i==7?'\n':' ');
    std::cout << std::endl;

    return 0;
}
