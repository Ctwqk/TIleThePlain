//
// Created by taiwei on 12/1/24.
//

#include "FourStarTwo.h"
int FourStarTwo::CHANNEL_=3;

FourStarTwo::FourStarTwo(Eigen::Vector3d conePoint1_, Eigen::Vector3d conePoint2_,int col, int row,std::string imgPath)
{
    if(imgPath.size()==0) imgPath="/home/taiwei/cg/Assignment3/code_part/data/1111.jpg";

    image = cv::imread(imgPath);
    int row_ = image.rows;
    int col_ = image.cols;
    if(row_>row||col_>col)
    {
        std::cerr<<"expecting size of view to be larger than size of tile"<<std::endl;
        return;
    }
    std::cout<<"start"<<std::endl;
    this->conePoint1 = conePoint1_; //left   (x,y,1)
    this->conePoint2 = conePoint2_; //right    (x,y,-1)
    this->conePoint1(0)+=(col-col_)/2;
    this->conePoint1(1)+=(row-row_)/2;
    this->conePoint2(0)+=(col-col_)/2;
    this->conePoint2(1)+=(row-row_)/2;
    this->row = row;
    this->col = col;
    angles[0]=M_PI/2;
    angles[1]=-M_PI/2;
    angles[2]=M_PI;
    angles[3]=0;

    // cv::imshow("image",image);
    // cv::waitKey(0);
    // cv::destroyAllWindows();
    tilePiece = new TilePiece(image,row,col);

    mediumTilePiece = new TilePiece(image,row,col);
    complement1=new TilePiece(cv::Mat::zeros(tilePiece->row,tilePiece->col,CV_8UC3),row,col); //left
    complement2=new TilePiece(cv::Mat::zeros(tilePiece->row,tilePiece->col,CV_8UC3),row,col); //right
    //std::cout<<complement1->Coos[0].size()<<std::endl;
    //std::vector<Eigen::MatrixXd> tmp(CHANNEL_,Eigen::MatrixXd::Zero(row,col)),tmp_(CHANNEL_,Eigen::MatrixXd::Zero(row,col));
    std::vector<Eigen::MatrixXd> tmp,tmp_,aux;

    //std::cout<<tmp[0].rows()<<" "<<tmp[0].cols()<<std::endl;
    int totalCol1 = 0, totalCol2 = 0;
    for(int i=0;i<3;i++)
    {
        //std::cout<<"rot"<<std::endl;
        aux=tilePiece->rotateCoo(conePoint1.head<2>(),angles[i]);
        tmp.insert(tmp.end(),aux.begin(),aux.end());
        totalCol1+=aux[0].cols();
        //complement1->Coos.insert(complement1->Coos.end(),tmp.begin(),tmp.end());

        //std::cout<<"rot_"<<std::endl;
        //std::cout<<"con"<<std::endl;
        //std::cout<<complement1->Coos[i].rows()<<" "<<complement1->Coos[i].cols()<<std::endl;
        //HelperFunctions::concate(tmp,complement1->imgs,tilePiece->row,tilePiece->col,0);

        aux=tilePiece->rotateCoo(conePoint2.head<2>(),angles[i]);
        tmp_.insert(tmp_.end(),aux.begin(),aux.end());
        totalCol2+=aux[0].cols();
        //complement2->Coos.insert(complement2->Coos.end(),tmp_.begin(),tmp_.end());
        //HelperFunctions::concate(tmp_,complement2->imgs,tilePiece->row,tilePiece->col,0);

        //std::cout<<"con_"<<std::endl;
    }
    complement1->Coos = std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(3,totalCol1));
    complement2->Coos = std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(3,totalCol1));
    //for(int i=0;i<CHANNEL_;i+=1)
    //{
        complement1->Coos[0]<<tmp[0],tmp[3],tmp[6];
        complement1->Coos[1]<<tmp[1],tmp[4],tmp[7];
        complement1->Coos[2]<<tmp[2],tmp[5],tmp[8];
        complement2->Coos[0]<<tmp_[0],tmp_[3],tmp_[6];
        complement2->Coos[1]<<tmp_[1],tmp_[4],tmp_[7];
        complement2->Coos[2]<<tmp_[2],tmp_[5],tmp_[8];
    //}
    // cv::Mat result;
    // HelperFunctions::cooToMat(complement1->Coos[0],tmp_[0],tilePiece->row,tilePiece->col);
    // HelperFunctions::cooToMat(complement1->Coos[2],tmp_[2],tilePiece->row,tilePiece->col);
    // HelperFunctions::cooToMat(complement1->Coos[1],tmp_[1],tilePiece->row,tilePiece->col);
    // HelperFunctions::eigenToMat(tmp_,result);
    //
    // cv::imshow("result",result);
    // cv::waitKey();
    // cv::destroyAllWindows();

    std::cout<<"end"<<std::endl;
}
int FourStarTwo::cols(){return col;}
int FourStarTwo::rows(){return row;}


void FourStarTwo::buildWithoutBfs()
{
    //TilePiece mediumTilePiece(cv::Mat::zeros(row,col,CV_8UC3));
    for(int i=0;i<CHANNEL_;i++)
    {

        mediumTilePiece->Coos[i]=Eigen::MatrixXd::Zero(complement2->Coos[i].rows(),complement2->Coos[i].cols()+tilePiece->Coos[i].cols());
        mediumTilePiece->Coos[i]<<tilePiece->Coos[i],complement2->Coos[i];
        //std::cout<<i<<std::endl;
    }
    //HelperFunctions::concate(mediumTilePiece->Coos,mediumTilePiece->imgs,row,col);
    std::vector<std::vector<Eigen::MatrixXd>> auxs(4);
    auxs[0] = tilePiece->rotateCoo(conePoint1.head<2>(),angles[0]);
    auxs[1] = tilePiece->rotateCoo(conePoint1.head<2>(),angles[1]);

    auxs[2] = tilePiece->rotateCoo(conePoint1.head<2>(),angles[2]);
    auxs[3]=std::vector<Eigen::MatrixXd>(3);
    HelperFunctions::transCoo(auxs[1],Eigen::Vector2d(400,0),auxs[3]);
    auxs[1]=auxs[3];
    HelperFunctions::transCoo(auxs[2],Eigen::Vector2d(200,-200),auxs[3]);
    auxs[2]=auxs[3];
    auxs[3] = tilePiece->transCoo(Eigen::Vector2d(200,200));
    HelperFunctions::concateCoo(auxs,mediumTilePiece->Coos);
    cv::Mat result;
    //for(int i=0;i<CHANNEL_;i++)HelperFunctions::cooToMat(mediumTilePiece->Coos[i],mediumTilePiece->imgs[i],row,col);
    HelperFunctions::concate(mediumTilePiece->Coos, mediumTilePiece->imgs,row,col);
    auxs[0]=mediumTilePiece->transCoo(Eigen::Vector2d(400,0));
    auxs[1]=mediumTilePiece->transCoo(Eigen::Vector2d(400,400));
    auxs[2]=mediumTilePiece->transCoo(Eigen::Vector2d(0,400));
    HelperFunctions::concate(auxs[0],mediumTilePiece->imgs,row,col);
    HelperFunctions::concate(auxs[1],mediumTilePiece->imgs,row,col);
    HelperFunctions::concate(auxs[2],mediumTilePiece->imgs,row,col);

    HelperFunctions::eigenToMat(mediumTilePiece->imgs,result);
    cv::imshow("result",result);
    cv::waitKey();
    cv::destroyAllWindows();

}

TilePiece* FourStarTwo::getMediumTilePiece(int LEVEL)
{

    int numInLayer;
    bool isOut = false;
    std::vector<Eigen::MatrixXd > tilesCoo =tilePiece->getCoos(), mTileCoo = mediumTilePiece->getCoos();
    Eigen::Vector3d cur,next;  //[2]>0: right, [2]<0: left
    Eigen::Vector2d disVec,disVec_,aux_;
    //double angles[]={M_PI/2,-M_PI/2,M_PI/2*3,0};
    std::vector<Eigen::Matrix2d> rotMats(4,Eigen::Matrix2d::Identity());
    rotMats[0] = HelperFunctions::getRotMatrix(angles[0],0);
    rotMats[1] = HelperFunctions::getRotMatrix(angles[1],0);
    rotMats[2] = HelperFunctions::getRotMatrix(angles[2],0);
    bfs1.push(conePoint1);
    bfs1.push(conePoint2);
    //Eigen::MatrixXd bfs1Coo(tilePiece->Coos[0]);
    std::vector<Eigen::MatrixXd> rotCoos(CHANNEL_),transCoos(CHANNEL_);


    aux_ = conePoint1.head<2>()- conePoint2.head<2>();
    std::vector<Eigen::Vector2d> points1={rotMats[0]*aux_+conePoint2.head<2>(),rotMats[1]*aux_+conePoint2.head<2>(),rotMats[2]*aux_+conePoint2.head<2>()};
    aux_ = -aux_;
    std::vector<Eigen::Vector2d> points2={rotMats[0]*aux_+conePoint1.head<2>(),rotMats[1]*aux_+conePoint1.head<2>(),rotMats[2]*aux_+conePoint1.head<2>()};
    std::vector<std::vector<Eigen::MatrixXd>> complements1={complement1->rotateCoo(conePoint1.head<2>(),angles[0],0,1),complement1->rotateCoo(conePoint1.head<2>(),angles[1],0,1),complement1->rotateCoo(conePoint1.head<2>(),angles[2],0,1),complement1->Coos};
    std::vector<std::vector<Eigen::MatrixXd>> complements2={complement2->rotateCoo(conePoint2.head<2>(),angles[0]),complement2->rotateCoo(conePoint2.head<2>(),angles[1]),complement2->rotateCoo(conePoint2.head<2>(),angles[2]),complement2->Coos};
    cv::Mat result;


    // HelperFunctions::cooToMat(complements1[3][0],rotCoos[0],tilePiece->row,tilePiece->col,0);
    // HelperFunctions::cooToMat(complements1[3][1],rotCoos[1],tilePiece->row,tilePiece->col,0);
    // HelperFunctions::cooToMat(complements1[3][2],rotCoos[2],tilePiece->row,tilePiece->col,0);
    // HelperFunctions::eigenToMat(rotCoos,result);
    // std::cout<<result.size()<<std::endl;
    // cv::imshow("result4",result);
    // cv::waitKey(0);
    // cv::destroyAllWindows();

    //std::cout<<rotCoos[0].rows()<<rotCoos[0].cols()<<std::endl;

    //std::cout<<complements1[0][0].rows()<<" "<<complements1[0][0].cols()<<std::endl;
    int angIdx,idx,layIdx=0;
    std::cout<<"start bfs"<<" "<<bfs1.size()<<std::endl;
    while(layIdx<LEVEL)
    {
        //isOut = true;
        numInLayer = bfs1.size();

        for(int k=0;k<numInLayer;k++)
        {
            cur = bfs1.front();
            bfs1.pop();
            if(cur[2]<0)
            {
                disVec = (cur-conePoint1).head<2>();
                //disVec_ = (cur-conePoint2).head<2>();
                angIdx = -cur[2]-1;
                for(int i=0;i<3;i++)
                {
                    idx = (angIdx+i)%3;
                    next<<points2[idx]+disVec,idx+1;
                    //std::cout<<next[0]<<" "<<next[1]<<std::endl;
                    if(next[0]<col&&next[0]>0&&next[1]<row&&next[1]>0)
                    {
                        bfs1.push(next);
                        //isOut = false;
                    }
                    else
                    {
                        bfs2.push(next);
                        bfs1.push(next);
                    }
                }
                //disVec = (cur - conePoint1).head<2>();
                //std::cout<<"rot"<<std::endl;
                HelperFunctions::transCoo(complements1[angIdx],disVec,transCoos);
                //std::cout<<"rot_"<<std::endl;
            }
            else
            {
                //disVec_ = (cur-conePoint1).head<2>();
                disVec = (cur-conePoint2).head<2>();
                angIdx = cur[2]-1;
                for(int i=0;i<3;i++)
                {
                    idx=(angIdx+i)%3;
                    next<<points1[idx]+disVec,-idx-1;
                    //std::cout<<next[0]<<" "<<next[1]<<std::endl;
                    if(next[0]<col&&next[0]>0&&next[1]<row&&next[1]>0) {
                        bfs1.push(next);
                        //isOut = false;
                    }
                    else
                    {
                        bfs2.push(next);
                        bfs1.push(next);
                    }
                }
                //std::cout<<"rot"<<std::endl;
                HelperFunctions::transCoo(complements2[angIdx],disVec,transCoos);
                //std::cout<<"rot_"<<std::endl;
            }
            //std::cout<<"trans_con"<<std::endl;

            HelperFunctions::concate(transCoos,mediumTilePiece->imgs,row,col,0);

            //std::cout<<"trans_con_"<<std::endl;
        }
        layIdx++;
        //std::cout<<bfs1.size()<<" "<<bfs2.size()<<std::endl;
    }
    std::cout<<"bfs end"<<std::endl;
    return mediumTilePiece;
}
int FourStarTwo::visualize(int col,int row)
{
    //for(int i=0;i<CHANNEL_;i++)HelperFunctions::cooToMat(mediumTilePiece->Coos[i], mediumTilePiece->imgs[i], row,col,0);
    cv::Mat result;
    HelperFunctions::eigenToMat(mediumTilePiece->imgs,result);
    //std::cout<<mediumTilePiece->imgs[1].maxCoeff()<<std::endl;
    //cv::resize(result,result,cv::Size(col,row));
    // cv::imshow("result",result);
    // cv::waitKey();
    // cv::destroyAllWindows();
    HelperFunctions::showImage(result,row,col);
    return 0;
}

FourStarTwo::~FourStarTwo()
{
    std::cout<<"Done"<<std::endl;
}
