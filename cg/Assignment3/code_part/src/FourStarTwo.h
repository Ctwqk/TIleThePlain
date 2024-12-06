//
// Created by taiwei on 12/1/24.
//

#ifndef FOURSTARTWO_H
#define FOURSTARTWO_H

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include "HelperFunctions.h"
#include <Eigen/Sparse>
#include <cmath>
#include "TilePiece.h"

class FourStarTwo {
    std::string imgPath ;
    static int CHANNEL_;
    Eigen::Vector3d conePoint1;
    Eigen::Vector3d conePoint2;
    int col, row;   //of final tile
    TilePiece *tilePiece;
    TilePiece *complement1, *complement2;
    cv::Mat image;
    double angles[4];
    TilePiece* mediumTilePiece;
    std::queue<Eigen::Vector3d> bfs1,bfs2;
public:
    FourStarTwo(Eigen::Vector3d conePoint1, Eigen::Vector3d conePoint2, int col, int row,std::string imgPath="");
    int cols();
    int rows();
    TilePiece* getMediumTilePiece(int LEVEL = 4);
    void buildWithoutBfs();
    int visualize(int col,int row);
    ~FourStarTwo();
};



#endif //FOURSTARTWO_H
