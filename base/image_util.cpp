//
// Created by rokeabbey on 2020/12/21.
//

#include "image_util.hpp"
#include <opencv4/opencv2/highgui.hpp>
#include <iostream>
#include <algorithm>

namespace rk {



    /**
     *
     * @tparam Tp 数组元素的类型
     * @param array 需要查找元素的数组, 必须先按照升序排列
     * @param target 目标元素
     * @param size 数组的大小
     * @param comparator 数组大小的比较器 第一个参数小于 第二个参数 则返回-1 反之返回1,相等返回0
     * @return 返回target在array中的下标， 如果target不在array中，则返回 -(元素需要插入的位置 + 1)
     */
    template<typename Tp>
    int binarySearch(const Tp *array, const Tp target, const int size, int (*comparator)(const Tp *, const Tp *)) {
        int l = 0, r = size - 1;
        int index;
        Tp v;
        while (l <= r) {
            index = (l + r) >> 1;
            v = array[index];
            if (comparator(&target, &v) < 0) {
                r = index - 1;
            } else if (comparator(&target, &v) > 0) {
                l = index + 1;
            } else {
                for (; index > 0; index--) {
                    v = array[index - 1];
                    if (comparator(&v, &target)) {
                        break;
                    }
                }
                return index;
            }
        }
        /**
        if (comparator(array + index, &target) < 0) {
            for (; index + 1 < size
                   && comparator(array + index, array + index + 1) == 0; index++) {}
        } else {
            for (; index - 1 > 0
                   && comparator(array + index, array + index - 1) == 0; index--) {}
        }
         */


        index = comparator(array + index, &target) < 0 ? index + 1 : index;
        return -index - 1;


    }

    int defaultIntComparator(const int *a, const int *b) {
        return *a - *b;
    }

    int getOptimizedSize(int n) {
        int index = binarySearch<int>(DEFAULT_OPTIMIZED_TABLE, n, DEFAULT_TABLE_LENGTH, defaultIntComparator);
        index = index < 0 ? -(index + 1) : index;
        return DEFAULT_OPTIMIZED_TABLE[index];
    }

    template<typename Tp>
    void arraySet(Tp *data, Tp value, int n) {
        Tp *p = data;
        for (int i = 0; i < n; i++, p++) {
            *p = value;
        }
    }



    namespace pdh {
        void main1() {
            int N = 1650;
            int padTempl[N];
            PaddingHelper<uchar> ph(10, 20, 30, 40, 0);
            ph.generatePaddingTemplate(padTempl, N);
            for (int i = 0; i < N; i++) {
                std::printf("%4d", padTempl[i]);
                if ((i + 1) & 15) {
                    std::printf(",");
                } else {
                    std::printf("\n");
                }
            }
            if (N & 15) {
                std::printf("\n");
            }
        }

        void main2() {
            PaddingHelper<uchar> ph(10, 20, 30, 40, 0);
            cv::Mat a(100, 110, CV_8UC3, cv::Scalar::all(100));
            cv::imshow("before", a);
            cv::waitKey(0);
            ph.setPaddingValue(a);
            cv::imshow("after", a);
            cv::waitKey(0);
        }

        void main3() {
            cv::Mat a(50, 70, CV_64FC3, cv::Scalar::all(100));
            PaddingHelper<double> ph(10, 5, 3, 7, 255);
            ph.setPaddingValue(a);
            printMat<cv::Vec3d>(a, true);
        }

        void main4() {
            for (int i = 0; i < rk::DEFAULT_TABLE_LENGTH; i++) {
                std::cout << rk::DEFAULT_OPTIMIZED_TABLE[i] << std::endl;
            }
//            std::cout << rk::DEFAULT_OPTIMIZED_TABLE[5] << std::endl;
        }

        void main0() {
            main4();
        }
    }
}
