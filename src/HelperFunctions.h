//
// Created by taiwei on 11/30/24.
//

#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H
#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include <cmath>

class HelperFunctions {
public:
    static Eigen::MatrixXd coordinateMat;
    static bool isInited;
    static double little;
    static int init(int row,int col);
    static int matToEigen(const cv::Mat &image, std::vector<Eigen::MatrixXd> &eigenImage );
    static int eigenToMat(const std::vector<Eigen::MatrixXd> &eigenImage, cv::Mat &image);
    static int rotateEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const Eigen::Vector2d &axis, Eigen::MatrixXd &eigenRotate, bool angleType = 0); //0: angle is 0 to 2pi 1:angle is -180 to 180
    static int transEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const double &distance,Eigen::MatrixXd &eigenTrans, bool angleType = 0, bool edgeOpe=0); //0:don't wrap edge 1: wrap edge
    static int reflection(const Eigen::MatrixXd &eigenImage, const double &angle, const double &distance,Eigen::MatrixXd &eigenTrans, bool angleType = 0, bool edgeOpe=0); //0:don't wrap edge 1: wrap edge

};



#endif //HELPERFUNCTIONS_H
