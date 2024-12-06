//
// Created by taiwei on 11/26/24.
//

#include <Eigen/Dense>
#include <opencv2/opencv.hpp>
#include "HelperFunctions.h"
#include "TilePiece.h"
#include "FourStarTwo.h"
#include "StarTripleThree.h"

int main(){
  // cv::Mat tile = cv::imread("/home/taiwei/cg/Assignment3/code_part/data/tile.png",cv::IMREAD_UNCHANGED);
  // if (tile.empty()) {
  //   std::cerr << "Error: Could not load image!" << std::endl;
  //   return -1;
  // }
  // std::vector<Eigen::MatrixXd> tile_mats;
  //
  // int row = 1200,col = 1200;

  // cv::Mat tile = cv::imread("/home/taiwei/cg/Assignment3/code_part/data/1111.jpg");
  // cv::imshow("tar",tile);
  //
  // TilePiece tp(tile);
  // std::vector<Eigen::MatrixXd> tileMat,tileCoo;
  // tp.Coos=tp.transCoo(Eigen::Vector2d(100,100));
  //
  // for(int i=0;i<3;i++)HelperFunctions::cooToMat(tp.Coos[i],tp.imgs[i],tile.rows,tile.cols);
  // cv::Mat result;
  // HelperFunctions::eigenToMat(tp.imgs,result);
  // cv::imshow("result1",result);
  // tp.Coos = tp.rotateCoo(Eigen::Vector2d(249,360),M_PI/2);
  // for(int i=0;i<3;i++)HelperFunctions::cooToMat(tp.Coos[i],tp.imgs[i],tile.rows,tile.cols);
  // HelperFunctions::eigenToMat(tp.imgs,result);
  // cv::imshow("result2",result);
  // cv::waitKey();
  // cv::destroyAllWindows();


  //
  // std::cout<<HelperFunctions::matToEigen(tile,tile_mats)<<std::endl;
  // int numOfChannels = tile_mats.size();
  // std::vector<Eigen::MatrixXd> paddedTileMats(numOfChannels,Eigen::MatrixXd::Zero(row,col));
  // int row_=tile_mats[0].rows();
  // int col_=tile_mats[0].cols();
  // for(int i=0;i<numOfChannels;i++)
  // {
  //   paddedTileMats[i].block(row/2-row_/2,col/2-col_/2,row_,col_)=tile_mats[i];
  // }
  //
  // HelperFunctions::init(col,row);
  // cv::Mat show;
  // Eigen::MatrixXd aux=paddedTileMats[0];
  // std::vector<Eigen::MatrixXd> auxs;
  // auxs.push_back(aux);
  // HelperFunctions::eigenToMat(auxs,show);
  // cv::imshow("show",show);

  //HelperFunctions::transEigen(paddedTileMats[0],45,900,aux,1,1);
  //Eigen::MatrixXd aux_ = Eigen::MatrixXd(aux);
  //HelperFunctions::transEigen(aux_,45,-900,aux,1,1);
  //HelperFunctions::transEigen(paddedTileMats[0],Eigen::Vector2d(-col_/2-col/2,-row_/2-row/2),aux,1);
  //auxs[0]=aux;s
  //HelperFunctions::rotateEigen(aux,-60,Eigen::Vector2d(col/2+col_/2,row/2+row_/2),aux_,1);
  //auxs[0]=aux_;
  //HelperFunctions::transEigen(aux_,Eigen::Vector2d(col/2+col_/2,row/2+row_/2),aux,1);
  //auxs[0]=aux;

  //FourStarTwo fourStarTwo(Eigen::Vector3d (149,260,-4), Eigen::Vector3d (350,260,4), 1000, 1000,"/home/taiwei/cg/Assignment3/code_part/data/1111.jpg");

  //for 4*2 and 442
  FourStarTwo fourStarTwo(Eigen::Vector3d (188,334,-4), Eigen::Vector3d (290,235,4), 600, 600,"/home/taiwei/cg/Assignment3/code_part/data/2222.jpg");
  fourStarTwo.getMediumTilePiece(5);
  fourStarTwo.buildWithoutBfs();
  fourStarTwo.visualize(800,800);

  //((222,317,0),(306,266,0),(224,214,0))
  //for *333
  // StarTripleThree starTripleThree(Eigen::Vector3d (222,317,0), Eigen::Vector3d (306,266,3), Eigen::Vector3d (224,214,6), 600, 600,"/home/taiwei/cg/Assignment3/code_part/data/3333.jpg");
  // starTripleThree.getMediumTilePiece(5);
  // starTripleThree.visualize(800,800);

  return 0;
}


