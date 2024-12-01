//
// Created by taiwei on 11/26/24.
//

#include "main.h"
#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include "HelperFunctions.h"

int main(){
  cv::Mat tile = cv::imread("/home/taiwei/cg/Assignment3/code_part/tile.png",cv::IMREAD_UNCHANGED);
  if (tile.empty()) {
    std::cerr << "Error: Could not load image!" << std::endl;
    return -1;
  }
  std::vector<Eigen::MatrixXd> tile_mats;

  int row = 2000,col = 2000;


  std::cout<<HelperFunctions::matToEigen(tile,tile_mats)<<std::endl;
  int numOfChannels = tile_mats.size();
  std::vector<Eigen::MatrixXd> paddedTileMats(numOfChannels,Eigen::MatrixXd::Zero(row,col));
  int row_=tile_mats[0].rows();
  int col_=tile_mats[0].cols();
  for(int i=0;i<numOfChannels;i++)
  {
    paddedTileMats[i].block(row/2-row_/2,col/2-col_/2,row_,col_)=tile_mats[i];
  }

  HelperFunctions::init(row,col);
  Eigen::MatrixXd aux;
  //HelperFunctions::transEigen(tile_mats[0],30,300,aux,1,1);
  HelperFunctions::rotateEigen(paddedTileMats[0],60,Eigen::Vector2d(0,0),aux,1);
  std::vector<Eigen::MatrixXd> auxs;
  auxs.push_back(aux);
  cv::Mat result;
  HelperFunctions::eigenToMat(auxs,result);

  std::cout<<HelperFunctions::eigenToMat(tile_mats,tile)<<std::endl;

  cv::imshow("result",result);
  cv::imshow("tile", tile);
  cv::waitKey();
  cv::destroyAllWindows();

  cv::imwrite( "/home/taiwei/cg/Assignment3/code_part/HW3_data/HW3_data/tile_processed.png",tile);
  return 0;
}


