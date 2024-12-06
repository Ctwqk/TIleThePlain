//
// Created by taiwei on 11/30/24.
//

#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H
#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <GLFW/glfw3.h>
#include <cmath>

class HelperFunctions {
public:
    static Eigen::MatrixXd coordinateMat;
    static bool isInited;
    static double little;
    static int CHANNEL_;
    static int init(int col,int row);
    static Eigen::Matrix2d getRotMatrix(double angle,bool angleType=0);
    static int matToEigen(const cv::Mat &image, std::vector<Eigen::MatrixXd> &eigenImages );
    static int eigenToMat(const std::vector<Eigen::MatrixXd> &eigenImages, cv::Mat &image);
    static int matToEigen(const cv::Mat &image, Eigen::MatrixXd &eigenImage );
    static int eigenToMat(const Eigen::MatrixXd &eigenImage, cv::Mat &image);
    static int rotateEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const Eigen::Vector2d &axis, Eigen::MatrixXd &eigenRotate, bool angleType = 0); //0: angle is 0 to 2pi 1:angle is -180 to 180
    static int transEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const double &distance,Eigen::MatrixXd &eigenTrans, bool angleType = 0, bool edgeOpe=0); //0:don't wrap edge 1: wrap edge
    static int transEigen(const Eigen::MatrixXd &eigenImage, Eigen::Vector2d disVec,Eigen::MatrixXd &eigenTrans, bool edgeOpe);
    static int reftEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const Eigen::Vector2d &startPoint, Eigen::MatrixXd &eigenRefted, bool angleType = 0, bool edgeOpe =0 );
    static int matToCoo(const Eigen::MatrixXd &eigenImage, Eigen::MatrixXd &Coo);
    static int cooToMat(const Eigen::MatrixXd &Coo, Eigen::MatrixXd &mat,int row,int col, bool edgeOpe= 0);
    static int normalize(const Eigen::MatrixXd &eigenImage, Eigen::MatrixXd &normalizedEigenImage);
    //static int padding(const Eigen::MatrixXd &eigenImage, Eigen::MatrixXd &paddedEigenImage,int row,int col,int padVal=0);
    static int concate(const std::vector<Eigen::MatrixXd> &Coos, std::vector<Eigen::MatrixXd> &resultEigens, int row,int col, bool edgeOpe=0);
    static int concateCoo( std::vector<std::vector<Eigen::MatrixXd>> &Coos,std::vector<Eigen::MatrixXd>& results);
    static int rotateCoo(std::vector<Eigen::MatrixXd>& imgCoos,  Eigen::Vector2d axis, double angle, std::vector<Eigen::MatrixXd> &rotateCoos, bool angleType=0 , bool edgeOpe=0 );
    static int transCoo(std::vector<Eigen::MatrixXd> &imgCoos, Eigen::Vector2d disVec,  std::vector<Eigen::MatrixXd> &transCoos,bool edgeOpe=0);
    static GLuint matToTex(const cv::Mat&mat);
    static int showImage(const cv::Mat&mat, int row, int col );
};



#endif //HELPERFUNCTIONS_H
