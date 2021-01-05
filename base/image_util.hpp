//
// Created by rokeabbey on 2020/12/21.
//

#ifndef FFT_IMPL_IMAGE_UTIL_HPP
#define FFT_IMPL_IMAGE_UTIL_HPP

#include <opencv4/opencv2/core.hpp>
#include <vector>
#include <assert.h>
#include <iostream>
//#include <algorithm>

namespace rk {


    template<typename PAD_TYPE>
    class PaddingHelper {

    public:
        PaddingHelper(int topBottom, int leftRight, PAD_TYPE paddingValue) : top(topBottom), bottom(topBottom),
                                                                             left(leftRight),
                                                                             right(leftRight),
                                                                             paddingValue(paddingValue) {}

        PaddingHelper(int top, int bottom, int left, int right, PAD_TYPE paddingValue) : top(top), bottom(bottom),
                                                                                         left(left), right(right),
                                                                                         paddingValue(paddingValue) {}

        void paddingImage(cv::Mat &src, cv::Mat &dst);

        static void generatePaddingTemplate(int padTempl[], int n);

        static int *generatePaddingTemplateWithNewArray(int n);
        static int *generatePaddingTemplateIncludeN(int n);
        void setPaddingValue(cv::Mat img);


    private:
        static int getMinAndChangeIndice(int pro[], int indice[]);

        int top, bottom, left, right;
        PAD_TYPE paddingValue;

    };

    template<typename Tp>
    void arraySet(Tp *data, Tp value, int n);

    static int DEFAULT_TABLE_LENGTH = 512;
    static int *DEFAULT_OPTIMIZED_TABLE = PaddingHelper<double>::generatePaddingTemplateWithNewArray(
            DEFAULT_TABLE_LENGTH);

    int getOptimizedSize(int n);

    template<typename Tp>
    int binarySearch(const Tp *array, const Tp target, const int size, int(*comparator)(const Tp *, const Tp *));

    int defaultIntComparator(const int *a, const int *b);

    template<typename Tp>
    void printMat(cv::Mat mat, bool printLineNum);

    template<typename PAD_TYPE>
    void PaddingHelper<PAD_TYPE>::setPaddingValue(cv::Mat img) {
        int nChannels = img.channels();
        int rowStep = nChannels * img.cols;
        int leftElemCount = nChannels * left;
        int rightElemCount = nChannels * right;
        int rightOffset = nChannels * (img.cols - right);
        bool bIsUCHAR = typeid(paddingValue).name() == typeid((uchar) 0).name();
        bool bUseMemset = bIsUCHAR;
        PAD_TYPE *data;
        data = (PAD_TYPE *) img.data;
        if (top > 0) {
            /** 最顶上的padding */
            PAD_TYPE *p = data;
            data += top * rowStep;
            if (bUseMemset) {
                memset(img.data, paddingValue, data - p);
            } else {
                arraySet<PAD_TYPE>(p, paddingValue, data - p);
            }
        }

        for (int r = top; r < img.rows - bottom; r++, data += rowStep) {
            /** 左边和右边的padding */
            if (bUseMemset) {
                memset(data, paddingValue, leftElemCount);
                memset(data + rightOffset, paddingValue, rightElemCount);
            } else {
                arraySet<PAD_TYPE>(data, paddingValue, leftElemCount);
                arraySet<PAD_TYPE>(data + rightOffset, paddingValue, rightElemCount);
            }
        }

        if (bottom > 0) {
            if (bUseMemset) {
                memset(data, paddingValue, bottom * rowStep);
            } else {
                arraySet<PAD_TYPE>(data, paddingValue, bottom * rowStep);
            }
        }
    }

    template<typename PAD_TYPE>
    void PaddingHelper<PAD_TYPE>::generatePaddingTemplate(int *padTempl, int n) {
        assert(n >= 1);
        padTempl[0] = 1;
        int indice[] = {0, 0, 0};
        for (int cur = 1; cur < n; cur++) {
            padTempl[cur] = getMinAndChangeIndice(padTempl, indice);
//            std::cout << "cur:" << cur << ", value:" << padTempl[cur] << std::endl;
        }
    }

    template<typename PAD_TYPE>
    int *PaddingHelper<PAD_TYPE>::generatePaddingTemplateWithNewArray(int n) {
        int *table = new int[n];
        generatePaddingTemplate(table, n);
        return table;
    }
    /**
     * 暂时不写了，核心代码更重要
     * @tparam PAD_TYPE
     * @param n
     * @return
     */
    template<typename PAD_TYPE>
    int *PaddingHelper<PAD_TYPE>::generatePaddingTemplateIncludeN(int n) {
        assert(n >= 1);
        int length = DEFAULT_TABLE_LENGTH;
        int *padTemp = new int[length];
        padTemp[0] = 1;
        int indice[] = {0, 0, 0};
        int cur = 1;
        do {
            for (; cur < length; cur ++){
                padTemp[cur] = getMinAndChangeIndice(padTemp, indice);
            }
        } while (n > padTemp[length - 1]);
        DEFAULT_OPTIMIZED_TABLE = padTemp;
        DEFAULT_TABLE_LENGTH = length;
    }

    template<typename PAD_TYPE>
    int PaddingHelper<PAD_TYPE>::getMinAndChangeIndice(int *padTempl, int *indice) {
        int pro2 = padTempl[indice[0]] << 1;
        int pro3 = padTempl[indice[1]] * 3;
        int pro5 = padTempl[indice[2]] * 5;
        int minValue = std::min(std::min(pro2, pro3), pro5);
        if (minValue == pro2) {
            indice[0]++;
        }

        if (minValue == pro3) {
            indice[1]++;
        }

        if (minValue == pro5) {
            indice[2]++;
        }
        return minValue;
    }

    /**
     *
     * @tparam PAD_TYPE
     * @param src
     * @param dst dst要与src不同哦
     */
    template<typename PAD_TYPE>
    void PaddingHelper<PAD_TYPE>::paddingImage(cv::Mat &src, cv::Mat &dst) {
        int dstColsCount = src.cols + left + right;
        int dstRowsCount = src.rows + top + bottom;
        dst.create(dstRowsCount, dstColsCount, src.type());
        dst = 0;
        cv::Mat a(10, 11, 12, 13);
        cv::Rect roiRect(left, top, src.cols, src.rows);
        src.copyTo(dst(roiRect));
        setPaddingValue(dst);
    }

    template<typename Tp>
    void printMat(cv::Mat mat, bool printLineNum) {
        for (int r = 0; r < mat.rows; r++) {
            Tp *p = mat.ptr<Tp>(r);
            for (int c = 0; c < mat.cols; c++, p++) {
                if (c) {
                    std::cout << ", ";
                } else if (printLineNum) {
                    std::cout << r << ": ";
                }
                std::cout << *p;
            }
            std::cout << std::endl;
        }
    }

    namespace pdh {
        void main0();
    }
}


#endif //FFT_IMPL_IMAGE_UTIL_HPP
