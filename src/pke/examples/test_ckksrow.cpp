// test_ckksrow.cpp
// Build: CMake / C++17 이상
// 목적: CKKS로 512차원 벡터 × (num_cols명) 내적을 행 기준 배치.
//       - 각 행(특징 i)의 전 사람 값을 CKKS 슬롯(= N/2)에 패킹(슬롯 초과 시 청크 분할)
//       - 질의 벡터 q[i]를 모든 슬롯에 복제한 평문과 ct-pt 곱: vec_len × chunks
//       - 곱 결과를 청크별로 누적(ct+ct) 후 복호
//       - FP32 데이터에 대해 표준화/정규화를 선택 적용 (행별 통계 사용)
//
// 실행 예:
//   (표준화; 기본) ./test_ckksrow 16384 512 8129 3 45 1 50
//   (정규화[-1,1]) ./test_ckksrow 16384 512 8129 3 45 1 50 2
// 인자: N vec cols runs scaleBits depth dcrtBits [preproc]
//   preproc: 0=none, 1=standardize(z-score, default), 2=normalize([-1,1])
//
// 주의:
//  - CKKS 슬롯 = N/2. num_cols > slots면 자동으로 여러 청크로 분할.
//  - scaleBits=45 권장(FP32 정밀도 여유), depth=1, dcrtBits≈50, FirstModSize≈60.

#include "openfhe.h"
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
#include <ios>
#include <iomanip>
#include <complex>
#include <cmath>
#include <limits>
#include <numeric>

using namespace lbcrypto;

template <class F>
static double time_ms(F&& fn) {
    auto t0 = std::chrono::high_resolution_clock::now();
    fn();
    auto t1 = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

static bool isPow2(size_t x){ return x && ((x & (x-1))==0); }

struct RowStats {
    double mean = 0.0;
    double std  = 1.0;
    double minv = 0.0;
    double maxv = 0.0;
};

enum class Preproc : int { None=0, Standardize=1, Normalize=2 };

static RowStats compute_row_stats(const std::vector<float>& rowF32) {
    RowStats st{};
    const size_t n = rowF32.size();
    if (n == 0) return st;

    // mean
    double sum = 0.0;
    double mn = std::numeric_limits<double>::infinity();
    double mx = -std::numeric_limits<double>::infinity();
    for (float xf : rowF32) {
        double x = static_cast<double>(xf);
        sum += x;
        if (x < mn) mn = x;
        if (x > mx) mx = x;
    }
    st.mean = sum / static_cast<double>(n);
    st.minv = mn;
    st.maxv = mx;

    // std (population std)
    double var = 0.0;
    for (float xf : rowF32) {
        double d = static_cast<double>(xf) - st.mean;
        var += d * d;
    }
    var /= static_cast<double>(n);
    st.std = std::sqrt(std::max(var, 0.0));
    return st;
}

static inline double apply_preproc(double x, const RowStats& st, Preproc mode) {
    constexpr double EPS = 1e-6;
    switch (mode) {
        case Preproc::Standardize: {
            double denom = std::max(st.std, EPS);
            return (x - st.mean) / denom;
        }
        case Preproc::Normalize: {
            double range = st.maxv - st.minv;
            double denom = std::max(range, EPS);
            return 2.0 * ((x - st.minv) / denom) - 1.0; // [-1,1]
        }
        case Preproc::None:
        default:
            return x; // no change
    }
}

int main(int argc, char** argv) {
    // args: N vec_len num_cols runs [scaleBits] [depth] [dcrtBits] [preproc]
    uint32_t N         = (argc > 1)? std::stoul(argv[1])  : 8192;   // ring dimension
    size_t   vec_len   = (argc > 2)? std::stoull(argv[2]) : 512;    // 차원(행 수)
    size_t   num_cols  = (argc > 3)? std::stoull(argv[3]) : 8129;   // 사람 수(열 수)
    int      runs      = (argc > 4)? std::stoi(argv[4])   : 1;      // 반복 횟수
    int      scaleBits = (argc > 5)? std::stoi(argv[5])   : 45;     // CKKS 스케일(log2 Δ)
    int      depth     = (argc > 6)? std::stoi(argv[6])   : 1;      // 곱 깊이
    int      dcrtBits  = (argc > 7)? std::stoi(argv[7])   : 50;     // 체인 소수 비트
    Preproc  preproc   = (argc > 8)? static_cast<Preproc>(std::stoi(argv[8])) : Preproc::Standardize;

    if (!isPow2(N)) {
        std::cerr << "[error] N must be a power of two.\n";
        return 1;
    }
    if (vec_len < 1) {
        std::cerr << "[error] vec_len must be ≥ 1.\n";
        return 1;
    }
    if (num_cols < 1) {
        std::cerr << "[error] num_cols must be ≥ 1.\n";
        return 1;
    }

    // ---- CKKS 컨텍스트 생성 (수동 파라미터) ----
    CCParams<CryptoContextCKKSRNS> ps;
    ps.SetSecurityLevel(HEStd_NotSet);  // 링 차원 직접 지정
    ps.SetRingDim(N);
    ps.SetMultiplicativeDepth(depth);   // ct*pt 1회
    ps.SetScalingModSize(dcrtBits);     // 체인의 (대부분) 소수 비트
    ps.SetFirstModSize(60);             // 초기 모듈러스(여유)

    auto cc = GenCryptoContext(ps);
    cc->Enable(PKE);
    cc->Enable(LEVELEDSHE);

    const size_t N_eff = cc->GetRingDimension();
    const size_t slots = N_eff / 2;                       // CKKS 슬롯 수(복소)
    const size_t chunks = (num_cols + slots - 1) / slots; // 청크 수

    // 모듈러스 체인 합 비트수(대략적인 log2 Q)
    auto ep = cc->GetCryptoParameters()->GetElementParams();
    using ParamsT = DCRTPoly::Params;
    auto dcrtParams = std::dynamic_pointer_cast<ParamsT>(ep);
    size_t numPrimes = dcrtParams ? dcrtParams->GetParams().size() : 0;
    double logQ = 0.0;
    if (dcrtParams){
        for (const auto& p : dcrtParams->GetParams())
            logQ += static_cast<double>(p->GetModulus().GetMSB());
    }

    const double scale = std::ldexp(1.0, scaleBits); // 2^scaleBits
    const char* preprocName = (preproc==Preproc::Standardize? "standardize(z-score)" :
                               preproc==Preproc::Normalize?   "normalize([-1,1])"   :
                                                             "none");

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "=== CKKS Row-Batched Dot (with FP32 preprocessing) ===\n"
              << "N="<<N<<" (eff "<<N_eff<<", slots="<<slots<<"), vec_len="<<vec_len
              << ", num_cols="<<num_cols<<", chunks="<<chunks<<", runs="<<runs<<"\n"
              << "scaleBits="<<scaleBits<<", depth="<<depth<<", dcrtBits="<<dcrtBits
              << ", preproc="<<preprocName
              << " -> numCRTprimes="<<numPrimes<<", approx log2(Q)="<<std::setprecision(1)<<logQ<<" bits\n\n";

    // ---- 키 생성 ----
    auto keys = cc->KeyGen();
    if (!keys.good()) {
        std::cerr << "KeyGen failed.\n";
        return 1;
    }

    // ---- 난수 생성기 (FP32 범위: 데모용으로 [-1000,1000]) ----
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<float> distFeat(-1000.0f, 1000.0f);
    std::uniform_real_distribution<float> distQuery(-1000.0f, 1000.0f);

    // ---- 행별 통계(μ/σ/min/max) 저장: vec_len 크기 ----
    std::vector<RowStats> stats(vec_len);

    // ---- DB 행 기준 배치(청크 분할) + 전처리 + 암호화 ----
    // ctRows[i][c] : i번째 행, c번째 청크의 암호문
    std::vector<std::vector<Ciphertext<DCRTPoly>>> ctRows(vec_len, std::vector<Ciphertext<DCRTPoly>>(chunks));

    for (size_t i=0;i<vec_len;++i) {
        // 1) 행 전체를 FP32로 생성 (num_cols 길이)
        std::vector<float> rowF32(num_cols);
        for (size_t j=0;j<num_cols;++j) rowF32[j] = distFeat(rng);

        // 2) 행 통계 계산
        stats[i] = compute_row_stats(rowF32);

        // 3) 청크 단위로 전처리 후 패킹/암호화
        for (size_t c=0;c<chunks;++c) {
            size_t base = c * slots;
            size_t used = std::min(slots, num_cols - base);

            std::vector<double> row(slots, 0.0); // CKKS 인코딩용(double)
            for (size_t j=0;j<used;++j) {
                double x = static_cast<double>(rowF32[base + j]);
                row[j] = apply_preproc(x, stats[i], preproc);
            }

            Plaintext ptRow = cc->MakeCKKSPackedPlaintext(row, scale);
            ctRows[i][c] = cc->Encrypt(keys.publicKey, ptRow);
        }
    }

    // ---- 워밍업 ----
    {
        // q[0] 생성 및 동일 행 통계(0번 행)로 전처리
        float q0f = distQuery(rng);
        double q0 = apply_preproc(static_cast<double>(q0f), stats[0], preproc);

        size_t c = 0;
        std::vector<double> s(slots, q0);
        Plaintext ptS = cc->MakeCKKSPackedPlaintext(s, scale);
        auto warmMul = cc->EvalMult(ctRows[0][c], ptS);
        cc->RescaleInPlace(warmMul);
        auto warmAcc = warmMul;
        if (vec_len > 1) warmAcc = cc->EvalAdd(warmAcc, warmMul);
        Plaintext warmOut;
        cc->Decrypt(keys.secretKey, warmAcc, &warmOut);
        (void)warmOut;
    }

    // ---- 본 측정 ----
    double sum_ms_mul = 0.0;   // vec_len*chunks ct*pt 총 시간
    double sum_ms_add = 0.0;   // (vec_len-1)*chunks ct+ct 총 시간
    double sum_ms_dec = 0.0;   // 청크 복호 총 시간
    double sum_ms_total = 0.0; // 전체 합

    for (int r=0; r<runs; ++r) {
        // 질의 벡터 q (FP32) 생성 후 행별 통계로 전처리
        std::vector<double> qProc(vec_len);
        for (size_t i=0;i<vec_len;++i) {
            float qi_f = distQuery(rng);
            qProc[i] = apply_preproc(static_cast<double>(qi_f), stats[i], preproc);
        }

        // partial[i][c] : i번째 행, c번째 청크의 부분곱
        std::vector<std::vector<Ciphertext<DCRTPoly>>> partial(vec_len, std::vector<Ciphertext<DCRTPoly>>(chunks));

        // ct*pt (행×청크)
        double ms_mul = time_ms([&]{
            for (size_t i=0;i<vec_len;++i) {
                std::vector<double> s(slots, qProc[i]); // q[i] 복제
                Plaintext ptS = cc->MakeCKKSPackedPlaintext(s, scale);
                for (size_t c=0;c<chunks;++c) {
                    partial[i][c] = cc->EvalMult(ctRows[i][c], ptS);
                    cc->RescaleInPlace(partial[i][c]); // 곱 1회 후 레벨/스케일 정렬
                }
            }
        });

        // 청크별 누적합
        std::vector<Ciphertext<DCRTPoly>> acc(chunks);
        double ms_add = time_ms([&]{
            for (size_t c=0;c<chunks;++c) {
                acc[c] = partial[0][c];
                for (size_t i=1;i<vec_len;++i) {
                    acc[c] = cc->EvalAdd(acc[c], partial[i][c]);
                }
            }
        });

        // 청크별 복호
        double ms_dec = 0.0;
        std::vector<Plaintext> out(chunks);
        ms_dec += time_ms([&]{
            for (size_t c=0;c<chunks;++c) {
                cc->Decrypt(keys.secretKey, acc[c], &out[c]);
            }
        });

        sum_ms_mul   += ms_mul;
        sum_ms_add   += ms_add;
        sum_ms_dec   += ms_dec;
        sum_ms_total += (ms_mul + ms_add + ms_dec);

        // (옵션) 결과 sanity check: 앞 8명(real part)만 출력
        if (r == runs-1) {
            std::cout << std::setprecision(6);
            std::cout << "[Sanity] first 8 persons (real parts): ";
            int printed = 0;
            for (size_t c=0;c<chunks && printed<8;++c) {
                auto vals = out[c]->GetCKKSPackedValue(); // vector<complex<double>>
                size_t base = c * slots;
                size_t used = std::min(slots, num_cols - base);
                for (size_t j=0;j<used && printed<8;++j) {
                    std::cout << vals[j].real() << (printed==7?'\n':' ');
                    ++printed;
                }
            }
        }
    }

    // ---- 결과 출력 ----
    const double runs_d = static_cast<double>(runs);
    const double ops_mul = runs_d * static_cast<double>(vec_len) * static_cast<double>( (chunks>0)?chunks:1 );
    const double ops_add = runs_d * static_cast<double>( (vec_len>0?vec_len-1:0) ) * static_cast<double>( (chunks>0)?chunks:1 );

    std::cout << std::setprecision(3);
    std::cout << "\n[Workload per run]\n"
              << "  ct*pt multiplies : " << (vec_len*chunks)  << " (row-wise × chunks)\n"
              << "  ct+ct adds       : " << ((vec_len>0?vec_len-1:0)*chunks) << "\n"
              << "  decrypts         : " << chunks << "\n\n";

    std::cout << "[Timing]\n"
              << "  PlainMult  total per run : " << (sum_ms_mul / runs_d) << " ms"
              << "  (avg per ct*pt: " << (sum_ms_mul / ops_mul) << " ms)\n"
              << "  CtAdd       total per run : " << (sum_ms_add / runs_d) << " ms"
              << "  (avg per add : " << (sum_ms_add / ops_add) << " ms)\n"
              << "  Decrypt     total per run : " << (sum_ms_dec / runs_d) << " ms\n"
              << "  TOTAL       per run       : " << (sum_ms_total / runs_d) << " ms\n";

    std::cout << std::endl;
    return 0;
}
