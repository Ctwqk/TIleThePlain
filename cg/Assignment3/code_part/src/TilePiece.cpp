//
// Created by taiwei on 12/1/24.
//

#include "TilePiece.h"

int TilePiece::CHANNEL_ = 3;
std::vector<Eigen::MatrixXd> TilePiece::transCoo( Eigen::Vector2d disVec,bool edgeOpe)
{
    Eigen::Vector3d addOn;
    std::vector<Eigen::MatrixXd> ans(3);
    addOn << disVec, 0;
    for(int i=0;i<3;i++)
    {
        ans[i]=(Coos[i].colwise()+addOn);
    }
    return ans;
}
std::vector<Eigen::MatrixXd> TilePiece::rotateCoo(  Eigen::Vector2d axis, double angle, bool angleType , bool edgeOpe )
{
    std::vector<Eigen::MatrixXd> auxs;
    std::vector<Eigen::MatrixXd> ans(3);
    int idx=0;
    //std::cout<<Coos.size()<<std::endl;
    //std::cout<<"trans"<<std::endl;
    HelperFunctions::transCoo(Coos, -axis, ans, edgeOpe);
    //std::cout<<"trans_"<<std::endl;
    Eigen::Matrix3d rotMat= Eigen::Matrix3d::Identity();
    rotMat.block(0,0,2,2)= HelperFunctions::getRotMatrix(angle, angleType);
    for(Eigen::MatrixXd &imgCoo : Coos)
    {
        auxs.push_back( rotMat*ans[idx]);
        idx++;
    }
    HelperFunctions::transCoo(auxs, axis, ans, edgeOpe);
    return ans;
}

TilePiece::TilePiece(cv::Mat mat,int row_,int col_)
{
    row = mat.rows;
    col = mat.cols;
    //std::cout<<"in"<<std::endl;
    Eigen::Vector3d disVec(0,0,0);
    HelperFunctions::matToEigen(mat,imgs);
    for(int i=0;i<CHANNEL_;i++)
    {
        Coos.push_back(Eigen::MatrixXd());
        HelperFunctions::matToCoo(imgs[i],Coos[i]);
        //std::cout<<imgs[i].cols()<<" "<<imgs[i].rows()<<std::endl;
    }
    if(row_>0)
    {
        disVec = Eigen::Vector3d((col_-col)/2,(row_-row)/2,0);
        std::vector<Eigen::MatrixXd> auxEigen=imgs;
        imgs = std::vector<Eigen::MatrixXd>(CHANNEL_,Eigen::MatrixXd::Zero(row_,col_));
        for(int i=0;i<CHANNEL_;i++)
        {

            imgs[i].block(disVec[1],disVec[0],row,col) = auxEigen[i];
            Coos[i].colwise() += disVec;
        }
        return;
    }

    HelperFunctions::matToEigen(mat,imgs);
    // cv::imshow("res",mat);
    // cv::waitKey(0);
    // cv::destroyAllWindows();
    //std::cout<<"out"<<std::endl;

    //std::cout<<Coos[0].cols()<<" "<<Coos[0].rows()<<std::endl;
}

int reftCoo(const Eigen::MatrixXd &imgCoo, Eigen::MatrixXd &reftCoo);
std::vector<Eigen::MatrixXd> TilePiece::getCoos()
{
    return Coos;
}
std::vector<Eigen::MatrixXd> TilePiece::getImgs()
{
    return imgs;
}