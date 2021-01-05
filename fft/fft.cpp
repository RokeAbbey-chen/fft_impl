//
// Created by rokeabbey on 2020/12/21.
//

#include "fft.hpp"
#include "../base/complex.hpp"
#include <iostream>
#include "../base/image_util.hpp"
#include <cmath>
#include <ctime>
#include <random>

namespace rk {
    /**
     *
     * @param n n要求大于1
     * @param factors
     * @return 质因子的数量（合并多个2为一个质因子，比如2*2*2则合并为8）, 如果n <= 1,则index返回0
     */
    int factorize(int n, int *factors) {
        int index = 0;
        factors[index] = ((n - 1) ^ n) + 1 >> 1;
        if (factors[index] > 1) {
            n /= factors[index];
            index++;
        }
        for (int f = 3; f * f <= n; f += 2) {
            for (; n / f * f == n; n /= f) {
                factors[index++] = f;
            }
        }
        if (n > 1) {
            factors[index++] = n;
        }
        if (factors[index - 1] != 1) {
            factors[index++] = 1;
        }
        return index;
    }

    /**
     * cooley-turkey算法 (蝴蝶算法)
     * @param info 要求将info.factors中的各个factors按照从大到小排列，并且将其中所有的2以其乘积的形式
     * 放在factors[0]的位置, 一个标准的factors 形如: factors = {8, 11, 5, 5, 3, 3}
     */
    void CT(FFT::Info &info) {
        int n = info.n;
        int nf = info.nf;
        int *factors = info.factors;
        int tailFactorsProduct[nf + 1];

        int digits[nf];
        for (int i = 0, j = n;
             i < nf;
             j /= factors[i], i++) {
            digits[i] = 0;
            tailFactorsProduct[i] = j;
        }
        std::cout << std::endl;
        tailFactorsProduct[nf] = 1;
        bool bIs2InFactors = (factors[0] & 1) == 0;
        info.iTab[0] = 0;
        if (bIs2InFactors) {
            for (int i = 0, inv = 0; i < factors[0]; i++) {
                info.iTab[i] = inv * tailFactorsProduct[1];
                for (int addNum = factors[0] >> 1;
                     addNum > 0 &&
                     ((inv ^= addNum) & addNum) == 0;
                     addNum >>= 1) {}
            }
        }

        if (bIs2InFactors && nf > 1 || nf > 0) {
            int startFIndex = bIs2InFactors;
            int baseIndex = bIs2InFactors ? factors[0] : 0;
            int calNumber = tailFactorsProduct[bIs2InFactors];
            for (int j = 1, v = 0; j < calNumber; j++) {
                int k = startFIndex;
                v += tailFactorsProduct[k + 1];
                digits[k]++;
                for (; digits[k] == factors[k]; k++) {
                    digits[k] = 0;
                    digits[k + 1]++;
                    v += tailFactorsProduct[k + 2] - tailFactorsProduct[k];
                }

                if (bIs2InFactors) {
                    for (int i = 0; i < factors[0]; i++) {
                        info.iTab[baseIndex + i] = info.iTab[i] + v;
                    }
                    baseIndex += factors[0];
                } else {
                    info.iTab[j] = v;
                }
            }
        }


    }

    template<typename Tp>
    void getWaveKTable(Complex<Tp> *_table, Complex<Tp> *_table2, int n, int cols) {
        if (n <= 0) { return; }
        Complex<Tp> (*table)[cols] = (Complex<Tp>(*)[cols]) _table;
        Complex<Tp> (*table2)[cols] = (Complex<Tp>(*)[cols]) _table2;
        Complex<Tp> wave(1, 0);
        Complex<Tp> delta;
        if (n == 2) {
            delta.re = -1, delta.im = 0;
        } else {
            double theta(-2 * M_PI / n);
            delta.re = cos(theta), delta.im = sin(theta);
        }
        Complex<Tp> one(1, 0);

        for (int i = 0; i < n; i++) {
            table[0][i] = one;
            table2[0][i] = one;
        }

        for (int i = 1; i < n; i++) {
            table[i][0] = one;
            table[i][1] = table[i - 1][1] * delta;

            table2[i][0] = one;
            table2[i][1] = table[i][1].conj();
            for (int j = 2; j < n; j++) {
                table[i][j] = table[i][j - 1] * table[i][1];
                table2[i][j] = table[i][j].conj();
            }
        }

    }

    template<typename Tp>
    void getWaveKMap(FFT::Info &info) {
        for (int index = 0; index < info.nf; index++) {
            if (index == 0 || info.factors[index] != info.factors[index - 1]) {
                int size = info.factors[index];
                if ((size & 1) == 0) {
                    size = 2;
                }
                if (!info.waveMap.count(size)) {
                    Complex<Tp> *wave = new Complex<Tp>[size * size];
                    Complex<Tp> *invWave = new Complex<Tp>[size * size];
                    getWaveKTable<Tp>(wave, invWave, size, size);
                    info.waveMap.insert(std::pair<int, void *>(size, wave));
                    info.waveMap.insert(std::pair<int, void *>(-size, invWave));
                }
            }
        }
    }

    /**
     * 对数据进行一次fft的计算
     * @tparam Tp
     * @param data 需要迭代计算的数据
     * @param size data的长度（等于一行或者一列的长度）
     * @param factor 当前所用到的质因子。是之前fftInit中分解质因数得到的一组数据中的一个，存储在rInfo或者cInfo中
     * @param baseRegion 上一次迭代时的计算域，与factor相乘(记为region)就代表第二层循环计算数据的多少
     * @param step 下标移动的步长， 当一行一行计算数据的时候是1, 当一列一列计算数据的时候 是列数
     * @param _waveKTable 提前计算好的wave表格
     * @param wKTableCols 表格的宽， 因为这里用到强制类型转换将一维指针转化成二维的形式，需要这个数据
     * @param sign 代表waveU迭代的时候是顺时针迭代还是逆时针迭代，-1代表顺时针（fft的方向）, 1代表逆时针(ifft的方向)
     */
    template<typename Tp>
    void calculateWithFactor(Complex<Tp> *data, int size, int factor, int baseRegion, int step,
                             Complex<Tp> *_waveKTable, int wKTableCols, int sign = -1) {

        int subRegion = baseRegion;
        int region = subRegion * factor;
        Complex<Tp> (*waveKTable)[wKTableCols] = (Complex<Tp> (*)[wKTableCols]) _waveKTable;
        /**
         * buffer用于存储 F[u + 0]，F[u + subRegion]， F[u + 2*subRegion], ... , F[u + (factor - 1) * subRegion]的中间计算结果
         */
        Complex<Tp> buffer[factor];
        Complex<Tp> zero(0, 0);
        double dThetaU = sign * M_PI * 2 / (factor * subRegion);
        Complex<Tp> dDeltaU(cos(dThetaU), sin(dThetaU));
        int forTimes = size / region;
        for (int t = 0; t < forTimes; t++, data += region * step) {
            Complex<Tp> deltaU(1, 0);
            /**
             * i就是fft计算中的自变量， 即F(i)中的自变量。
             * fac = factor, sub = subRegion, region = factor * subRegion
             *
             * F(u) = F0(u) * 1 + F1(u) * exp{-j2pi * 1(u)/region} + F2(u) * exp{-j2pi * 2(u) / region} ....
             * F(u + sub) = F0(u + sub) * 1 + F1(u + sub) * exp{-j2pi * 1(u + sub) / region} + F2(u + sub) * exp{-j2pi * 2(u + sub)/region} ...
             *            = F0(u + sub)     + F1(u + sub) * exp{-j2pi·u/region} * exp{-j2pi·1/factor}
             *              + F2(u + 2sub) * exp{-j2pi·2u/region} * exp{-j2pi·2/factor}
             *              这里注意一下等式的变换，拆分出来了方便循环计算, 还有 region = factor * sub要注意
             * ...
             * ...
             * F(u + (fac - 1)sub) = F0(u + (fac - 1)sub) * exp{-j2pi * 1(u + (fac - 1)sub) / region}
             *                      + F2(u + (fac - 1)sub) * exp{-j2pi * 2(u + (fac - 1)sub)/region} ...
             *
             * i这层循环沿着 u in [0, subRegion)移动
             * k这层循环 沿着 F0, F1, F2...方向移动计算
             * j这层循环沿着 u, u + sub, u + 2sub这个方向移动计算
             *
             * 值得注意的是计算的时候有很多 复数共轭的性质可以使用，这里是因为先写好了代码，再去opencv中参考才发现可以使用共轭的性质简化计算，故没有使用
             *
             * 由于操作符重载或者调用函数计算都会有额外开辟栈帧的消耗，所以下面的计算全部换成宏函数，提高效率。
             */
            for (int i = 0; i < subRegion; i++) {
                for (int j = 0; j < factor; j++) {
                    buffer[j] = zero;
                }

                Complex<Tp> waveU(1, 0);
                for (int k = 0, offset = 0; k < factor; k++, offset += subRegion) {
                    int dataIndex = (i + offset) * step;
                    Complex<Tp> dataWaveU;//, tmp;
                    COMPLEX_MULTI(dataWaveU, data[dataIndex], waveU)

                    for (int j = 0; j < factor; j++) {
                        COMPLEX_ADD_ASSIGN_MULTI(buffer[j], dataWaveU, waveKTable[k][j])
                    }

                    Tp waveUTmp;
                    COMPLEX_MULTI_ASSIGN(waveU, deltaU, waveUTmp)
                }
                for (int j = 0, offset = 0; j < factor; j++, offset += subRegion) {
                    data[(i + offset) * step] = buffer[j];
                }
                Tp deltaUTmp;
                COMPLEX_MULTI_ASSIGN(deltaU, dDeltaU, deltaUTmp)
            }
        }
    }

    /**
     * 根据src的行数列数分解质因数的结果对每个质因数调用 calculateWithFactor 进行计算，
     * @tparam Tp
     * @param src 源矩阵起始指针, 根据此指针以及info中的size，可以对一行或者一列进行fft的计算
     * @param dst 目标矩阵起始指针, 根据次指针以及info中的size，可以对一行或者一列进行fft结果的存放
     * @param info
     */
    template<typename Tp>
    void doCalculate(Complex<Tp> *src, Complex<Tp> *dst, FFT::Info &info, int notInv = 1) {
        Complex<Tp> *_dst = dst;
        int dstStep = info.step;
        if (src == _dst) {
            dst = new Complex<Tp>[info.n];
            dstStep = 1;
        }

        for (int i = 0, index = 0; i < info.n; i++, index += dstStep) {
            dst[index] = src[info.iTab[i] * info.step]; // 0 2 1 3
        }
        int baseRegion = 1;
        Complex<Tp> *waveKTable;
        for (int i = 0; i < info.nf; i++) {
            int fct = info.factors[i];
            if ((fct & 1) == 0) {
                waveKTable = (Complex<Tp> *) info.waveMap.find(notInv * 2)->second;
                while ((fct & 1) == 0) {
                    calculateWithFactor<Tp>(dst, info.n, 2, baseRegion, dstStep, waveKTable, 2, -notInv);
                    baseRegion <<= 1;
                    fct >>= 1;
                }
            } else {
                if (i == 0 || info.factors[i] != info.factors[i - 1]) {
                    waveKTable = (Complex<Tp> *) info.waveMap.find(notInv * info.factors[i])->second;
                }
                calculateWithFactor<Tp>(dst, info.n, fct, baseRegion, dstStep, waveKTable, fct, -notInv);
                baseRegion *= fct;
            }
        }
        if (src == _dst) {
            for (int i = 0, index = 0; i < info.n; i++, index += info.step) {
                _dst[index] = dst[i];
            }
            delete[]dst;
        }
    }

    void FFT::fftInit() {
        dst.create(src.size(), src.type());

        rInfo.n = src.cols;
        cInfo.n = src.rows;

        rInfo.nf = factorize(src.cols, rInfo.factors);
        cInfo.nf = factorize(src.rows, cInfo.factors);

        rInfo.step = 1;
        cInfo.step = src.cols;


        rInfo.iTab = new int[rInfo.n];
        cInfo.iTab = new int[cInfo.n];

        if (src.depth() == CV_64F) {
//            rInfo.wave = getWave<double>()
        }

        CT(rInfo);
        CT(cInfo);

        /* 计算waveK*/
        if (src.depth() == CV_64F) {
            getWaveKMap<double>(rInfo);
            cInfo.waveMap = rInfo.waveMap;
            getWaveKMap<double>(cInfo);
        } else {
            getWaveKMap<float>(rInfo);
            cInfo.waveMap = rInfo.waveMap;
            getWaveKMap<float>(cInfo);
        }

    }


    template<typename Tp>
    cv::Mat &FFT::calc(int notInv, bool direct) {
        Info infos[] = {cInfo, rInfo};

        clock_t start = clock();
        for (int d = 0; d < 2; d++, direct = !direct) {
            Complex<Tp> *srcP = (Complex<Tp> *) src.data;
            Complex<Tp> *dstP = (Complex<Tp> *) dst.data;
            if (d == 1) {
                srcP = dstP;
            }
            Info curInfo = infos[direct];
            for (int i = 0; i < infos[!direct].n;
                 i++, srcP += infos[!direct].step, dstP += infos[!direct].step) {
                doCalculate<Tp>(srcP, dstP, curInfo, notInv);
            }
        }
        clock_t end = clock();
        std::cout << "核心计算耗时:" << (end - start) / 1000 << "ms\n";
        Complex<Tp> *dstP = (Complex<Tp> *) dst.data;
        for (int i = 0; i < dst.rows; i++) {
            for (int j = 0; j < dst.cols; j++, dstP++) {
                if (dstP->re < 1e-10 && dstP->re > -1e-10) {
                    dstP->re = 0;
                }
                if (dstP->im < 1e-10 && dstP->im > -1e-10) {
                    dstP->im = 0;
                }
            }
        }

        return dst;
    }

    cv::Mat &FFT::calc(int notInv, bool direct) {
        if (src.depth() == CV_64F) {
            return calc<double>(notInv, direct);
        } else if (src.depth() == CV_32F) {
            return calc<float>(notInv, direct);
        } else {
            throw "请使用double / float类型的数据 （CV_64F / CV_32F）";
        }
    }

    namespace fft_demo {
        void main8() {
            int rows = 3 * 5;// 4 * 9 * 12 * 5;
            int cols = 3 * 5;// 4 * 9 * 12 * 5;
            cv::Mat src = cv::Mat::zeros(rows, cols, CV_64FC2);
            srand(time(0));
            std::cout << "rand() : " << rand() << std::endl;
            for (int k = 0; k < rows; k++) {
                for (int i = 0; i < cols; i++) {
                    src.at<cv::Vec2d>(k, i) = cv::Vec2d(k * cols + (i + 1) + (rand() & 15), 0);
                }
            }

            std::cout << "---------原Mat输出--------------" << std::endl;
            printMat<cv::Vec2d>(src, false);

            cv::Mat dst, dst2;
            dst2.create(cv::Size(10, 10), src.depth());
            FFT fft(src, dst);
            fft.fftInit();
            fft.calc(1, true);
            std::cout << "--------- 本算法fft 结果输出 --------------" << std::endl;
            printMat<cv::Vec2d>(dst, false);

            FFT ifft(dst, dst2);
            ifft.fftInit();
            ifft.calc(-1, true);
            std::cout << "--------- 本算法 ifft 结果输出--------------" << std::endl;
            printMat<cv::Vec2d>(dst2, false);


            cv::Mat dst3;
            cv::idft(dst, dst3);
            std::cout << "--------- opencv 4.5.0 ifft输出--------------" << std::endl;
            printMat<cv::Vec2d>(dst3, false);

        }

        void main0() {
            main8();
        };


    }
}
