//
// Created by rokeabbey on 2020/12/21.
//

#ifndef FFT_IMPL_FFT_HPP
#define FFT_IMPL_FFT_HPP

#include<opencv4/opencv2/core.hpp>

#include<map>

namespace rk {

    class FFT {

    public:
        class Info {
        public:
            /**
             * 用于存储图像的行长或者列长, 也是factors所有因子的乘积
             */
            int n;
            int factors[34];
            /**
             * length of the factors
             */
            int nf;
            /**
             * 元素的计算间隔，当沿着列方向计算的时候为1,当沿着行方向计算的时候为列数
             */
            int step;

            int *iTab = NULL;

            std::map<int, void *> waveMap;


        };

    private:
        bool needPadding;
        cv::Mat &src, &dst;
        Info rInfo, cInfo;


    public:
        FFT(cv::Mat &src, cv::Mat &dst) : src(src), dst(dst), needPadding(true) {}

        /**
         *
         * @param src 输入图像，需要是复数形式，即矩阵包含两个通道，第一个通道为图像像素的实部，第二个通道为像素的虚部（一般直接读进来的图像虚部设为0）
         * @param dst 输出图像，也是两个通道，一个实部，一个虚部
         * @param needPadding
         */
        FFT(cv::Mat &src, cv::Mat &dst, bool needPadding) : src(src), dst(dst), needPadding(needPadding) {}


        void fftInit();

        /**
         * @param notInv 1 表示计算fft， -1表示计算ifft.
         * @param direct true代表沿着列增长方向计算，false代表沿着行方向计算
         * */
        cv::Mat &calc(int notInv = 1, bool direct = true);

        template<typename name>
        cv::Mat &calc(int notInv = 1, bool direct = true);


    };
    namespace fft_demo {
        void main0();
    }

}


#endif //FFT_IMPL_FFT_HPP
