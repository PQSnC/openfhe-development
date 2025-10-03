//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//==================================================================================

/*
  Example of Proxy Re-Encryption on a packed vector.
  Multiparty proxy-reencryption of an integer buffer using BFV RNS scheme,
  now with command-line tunables.
*/

#define PROFILE  // for TIC/TOC timing
#include "openfhe.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace lbcrypto;

using CT = Ciphertext<DCRTPoly>;
using PT = Plaintext;

using vecInt  = std::vector<int64_t>;
using vecChar = std::vector<char>;

static bool run_demo_pre(uint32_t N,
                         size_t vec_len,
                         size_t num_cols,
                         int runs,
                         uint64_t t,
                         int depth,
                         int dcrtBits,
                         int ksArg,
                         int mulArg,
                         int measureRot,
                         int autoMode,
                         long long kConst);

//----------------------------- MAIN -------------------------------------------

int main(int argc, char** argv) {
    // args: N vec_len num_cols runs [t] [depth] [dcrtBits] [ks] [mul] [measureRot] [autoMode] [kConst]
    uint32_t   N          = (argc > 1)  ? static_cast<uint32_t>(std::stoul(argv[1]))  : 2048;
    size_t     vec_len    = (argc > 2)  ? static_cast<size_t>(std::stoull(argv[2]))   : 512;
    size_t     num_cols   = (argc > 3)  ? static_cast<size_t>(std::stoull(argv[3]))   : 2048;
    int        runs       = (argc > 4)  ? std::stoi(argv[4])                           : 1;
    uint64_t   t          = (argc > 5)  ? std::stoull(argv[5])                         : 65537;
    int        depth      = (argc > 6)  ? std::stoi(argv[6])                           : 1;
    int        dcrtBits   = (argc > 7)  ? std::stoi(argv[7])                           : 60;
    int        ksArg      = (argc > 8)  ? std::stoi(argv[8])                           : 0;   // 0:BV, 1:HYBRID, 2:GHS
    int        mulArg     = (argc > 9)  ? std::stoi(argv[9])                           : 0;   // 1이면 EvalMultKeyGen
    int        measureRot = (argc > 10) ? std::stoi(argv[10])                          : 0;   // 1이면 회전키 생성
    int        autoMode   = (argc > 11) ? std::stoi(argv[11])                          : 0;   // 1이면 N/depth 자동화
    long long  kConst     = (argc > 12) ? std::stoll(argv[12])                         : 5;

    bool ok = run_demo_pre(N, vec_len, num_cols, runs, t, depth, dcrtBits,
                           ksArg, mulArg, measureRot, autoMode, kConst);
    return ok ? 0 : 1;
}

//--------------------------- IMPLEMENTATION -----------------------------------

static bool run_demo_pre(uint32_t N,
                         size_t vec_len,
                         size_t num_cols,
                         int runs,
                         uint64_t t,
                         int depth,
                         int dcrtBits,
                         int ksArg,
                         int mulArg,
                         int measureRot,
                         int autoMode,
                         long long kConst) {
    // 타이머
    TimeVar timer;

    std::cout << "================ PRE(BFV-RNS) Demo =================\n";
    std::cout << "Params from CLI:\n"
              << "  N=" << N << "  vec_len=" << vec_len << "  num_cols=" << num_cols
              << "  runs=" << runs << "  t=" << t << "  depth=" << depth
              << "  dcrtBits=" << dcrtBits << "  ks=" << ksArg
              << "  mul=" << mulArg << "  measureRot=" << measureRot
              << "  autoMode=" << autoMode << "  kConst(seed)=" << kConst << "\n\n";

    // -------------------- CryptoContext 설정 --------------------
    std::cout << "Setting up BFV RNS crypto system..." << std::endl;

    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetPlaintextModulus(t);
    parameters.SetScalingModSize(dcrtBits);

    // 키 스위칭 기법 선택
    switch (ksArg) {
    case 1:
        parameters.SetKeySwitchTechnique(KeySwitchTechnique::HYBRID);
        break;
    case 0:
    default:
        // ksArg가 0이거나 알 수 없으면 BV로
        parameters.SetKeySwitchTechnique(KeySwitchTechnique::BV);
        if (ksArg != 0 && ksArg != 1) {
            std::cerr << "[WARN] ksArg=" << ksArg
                      << " is not supported. Using BV.\n";
        }
        break;
}

    if (autoMode) {
        // 자동 모드: 보안수준만 지정하고 나머지는 내부에서 결정
        parameters.SetSecurityLevel(SecurityLevel::HEStd_128_classic);
    } else {
        // 수동 모드: 링 차원과 곱셈 깊이 직접 지정
        parameters.SetRingDim(N);
        if (depth > 0) {
            parameters.SetMultiplicativeDepth(depth);
        }
    }

    TIC(timer);
    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    std::cout << "\nParam generation time:\t" << TOC_MS(timer) << " ms" << std::endl;

    // 기능 활성화
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(PRE);

    // 요약 정보 출력
    std::cout << "p (pt modulus) = " << cc->GetCryptoParameters()->GetPlaintextModulus() << std::endl;
    std::cout << "n (ring dim)   = "
              << cc->GetCryptoParameters()->GetElementParams()->GetCyclotomicOrder() / 2 << std::endl;
    std::cout << "log2(q)        = "
              << log2(cc->GetCryptoParameters()->GetElementParams()->GetModulus().ConvertToDouble()) << std::endl;
    std::cout << "r (digit size) = " << cc->GetCryptoParameters()->GetDigitSize() << std::endl;

    auto ringsize = cc->GetRingDimension();
    std::cout << "Slots available: " << ringsize << " (Alice can encrypt ~" << ringsize * 2 << " bytes)\n";

    if (vec_len > ringsize) {
        std::cout << "[WARN] vec_len(" << vec_len << ") > slots(" << ringsize
                  << "). vec_len을 " << ringsize << "로 클램프합니다.\n";
        vec_len = ringsize;
    }

    if (num_cols != 0) {
        std::cout << "[INFO] num_cols 파라미터는 본 PRE 데모에서는 사용하지 않습니다.\n";
    }

    // 난수 시드(재현성)
    std::mt19937_64 rng(static_cast<uint64_t>(kConst));
    std::uniform_int_distribution<int64_t> dist(0, static_cast<int64_t>(t - 1));

    bool all_good = true;

    for (int run = 1; run <= runs; ++run) {
        std::cout << "\n---- Run " << run << " / " << runs << " ----\n";

        // -------------------- Alice 키 생성 --------------------
        std::cout << "Alice key generation..." << std::endl;
        TIC(timer);
        KeyPair<DCRTPoly> keyPair1 = cc->KeyGen();
        std::cout << "KeyGen time:\t\t" << TOC_MS(timer) << " ms" << std::endl;
        if (!keyPair1.good()) {
            std::cout << "Alice Key generation failed!\n";
            return false;
        }

        if (mulArg) {
            // 본 데모에선 곱셈을 하진 않지만 옵션으로 키만 생성
            cc->EvalMultKeyGen(keyPair1.secretKey);
            std::cout << "EvalMultKeyGen done.\n";
        }

        if (measureRot) {
            // 회전 키 생성(예: ±1, ±2, ±4)
            std::vector<int32_t> idx{1, 2, 4, -1, -2, -4};
            cc->EvalAtIndexKeyGen(keyPair1.secretKey, idx);
            std::cout << "EvalAtIndexKeyGen for indices {±1,±2,±4} done.\n";
        }

        // -------------------- 데이터 인코딩 --------------------
        vecInt vShorts(vec_len);
        for (size_t i = 0; i < vec_len; ++i)
            vShorts[i] = dist(rng);

        PT pt = cc->MakePackedPlaintext(vShorts);

        // -------------------- 암호화 --------------------
        TIC(timer);
        auto ct1 = cc->Encrypt(keyPair1.publicKey, pt);
        std::cout << "Encrypt time:\t\t" << TOC_MS(timer) << " ms" << std::endl;

        // -------------------- 복호화(검증용) --------------------
        PT ptDec1;
        TIC(timer);
        cc->Decrypt(keyPair1.secretKey, ct1, &ptDec1);
        std::cout << "Decrypt time(Alice):\t" << TOC_MS(timer) << " ms" << std::endl;
        ptDec1->SetLength(pt->GetLength());

        // -------------------- Bob 키 생성 --------------------
        std::cout << "Bob key generation..." << std::endl;
        TIC(timer);
        KeyPair<DCRTPoly> keyPair2 = cc->KeyGen();
        std::cout << "KeyGen time:\t\t" << TOC_MS(timer) << " ms" << std::endl;
        if (!keyPair2.good()) {
            std::cout << "Bob Key generation failed!\n";
            return false;
        }

        // -------------------- 재암호화 키 생성 --------------------
        std::cout << "Generating re-encryption key..." << std::endl;
        TIC(timer);
        EvalKey<DCRTPoly> reencryptionKey12 = cc->ReKeyGen(keyPair1.secretKey, keyPair2.publicKey);
        std::cout << "ReKeyGen time:\t\t" << TOC_MS(timer) << " ms" << std::endl;

        // -------------------- 재암호화 --------------------
        TIC(timer);
        auto ct2 = cc->ReEncrypt(ct1, reencryptionKey12);
        std::cout << "Re-Encrypt time:\t" << TOC_MS(timer) << " ms" << std::endl;

        // -------------------- Bob 복호화 --------------------
        PT ptDec2;
        TIC(timer);
        cc->Decrypt(keyPair2.secretKey, ct2, &ptDec2);
        std::cout << "Decrypt time(Bob):\t" << TOC_MS(timer) << " ms" << std::endl;
        ptDec2->SetLength(pt->GetLength());

        // -------------------- 결과 검증 --------------------
        auto unpacked0 = pt->GetPackedValue();
        auto unpacked1 = ptDec1->GetPackedValue();
        auto unpacked2 = ptDec2->GetPackedValue();

        // OpenFHE는 평문을 [-p/2, p/2)로 해석하므로 음수면 p를 더해 0..p-1로 복원
        for (size_t j = 0; j < pt->GetLength(); ++j) {
            if (unpacked1[j] < 0) unpacked1[j] += static_cast<int64_t>(t);
            if (unpacked2[j] < 0) unpacked2[j] += static_cast<int64_t>(t);
        }

        bool good = true;
        for (size_t j = 0; j < pt->GetLength(); ++j) {
            if ((unpacked0[j] != unpacked1[j]) || (unpacked0[j] != unpacked2[j])) {
                std::cout << "Mismatch at " << j
                          << ": src=" << unpacked0[j]
                          << ", decA=" << unpacked1[j]
                          << ", decB=" << unpacked2[j] << std::endl;
                good = false;
                break;
            }
        }

        if (good) std::cout << "[Run " << run << "] PRE passes\n";
        else      std::cout << "[Run " << run << "] PRE fails\n";

        all_good = all_good && good;
    }

    std::cout << "\n================ Execution Completed ================\n";
    return all_good;
}
