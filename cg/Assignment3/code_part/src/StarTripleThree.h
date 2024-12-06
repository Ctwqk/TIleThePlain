//
// Created by taiwei on 12/1/24.
//

#ifndef STARTRIPLETHREE_H
#define STARTRIPLETHREE_H

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include "HelperFunctions.h"
#include <Eigen/Sparse>
#include <cmath>
#include "TilePiece.h"

class StarTripleThree {
    std::string imgPath ;
    static int CHANNEL_;
    Eigen::Vector3d conePoint1;
    Eigen::Vector3d conePoint2;
    Eigen::Vector3d conePoint3;
    int col, row;   //of final tile
    TilePiece *tilePiece;
    TilePiece *complement1, *complement2, *complement3;
    cv::Mat image;
    double angles[3];
    TilePiece* mediumTilePiece;
    std::queue<Eigen::Vector3d> bfs1,bfs2;
public:
    StarTripleThree(Eigen::Vector3d conePoint1, Eigen::Vector3d conePoint2,Eigen::Vector3d conePoint3,int col, int row,std::string imgPath="");
    int cols();
    int rows();
    TilePiece* getMediumTilePiece(int LEVEL = 4);
    int visualize(int col,int row);
    ~StarTripleThree();
};



#endif //FOURSTARTWO_H
