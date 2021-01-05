//
// Created by rokeabbey on 2020/12/22.
//

#include "boot.hpp"
#include "opencv4/opencv2/core.hpp"
//#include "../fft/fft.hpp"
#include "../base/image_util.hpp"

#include <algorithm>

namespace rk {
    int comprator(const void *a, const void *b) {
        return *((int *) a) - *((int *) b);

    }

    namespace entrance {
        void main1() {
            int nRow = 101, nCol = 33;
            cv::Mat mat = cv::Mat::zeros(nRow, nCol, CV_64FC1);
            int optNRow = rk::getOptimizedSize(mat.rows);
            int optNCol = rk::getOptimizedSize(mat.cols);
            cv::Mat padded;
            rk::PaddingHelper<double> paddingHelper(0, optNRow - mat.rows, 0, optNCol - mat.cols, 1);
            paddingHelper.paddingImage(mat, padded);
            rk::printMat<cv::Vec<double, 1>>(padded, true);
        }

        void main2(){
            std::vector<double> a;
            a.pop_back();
        }

        void main0() {
            main1();
        }

    }
}