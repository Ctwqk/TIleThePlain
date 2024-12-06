//
// Created by taiwei on 12/1/24.
//

#ifndef TILEPIECE_H
#define TILEPIECE_H

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include "HelperFunctions.h"
#include <Eigen/Sparse>
#include <cmath>

class TilePiece {
public:
      static int CHANNEL_;
      std::vector<Eigen::MatrixXd> Coos;
      std::vector<Eigen::MatrixXd> imgs;
      int row;
      int col;
      std::vector<Eigen::MatrixXd> transCoo( Eigen::Vector2d disVec,  bool edgeOpe = 0);
      std::vector<Eigen::MatrixXd> rotateCoo( Eigen::Vector2d axis, double angle, bool angleType = 0, bool edgeOpe = 0 );
      int reftCoo(const Eigen::MatrixXd &imgCoo, Eigen::MatrixXd &reftCoo);
      TilePiece(cv::Mat mat,int row = -1,int col = -1);
      std::vector<Eigen::MatrixXd> getCoos();
      std::vector<Eigen::MatrixXd> getImgs();
};



#endif //TILEPIECE_H
