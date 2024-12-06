//
// Created by taiwei on 12/5/24.
//

#include "StarTripleThree.h"
int StarTripleThree::CHANNEL_=3;

StarTripleThree::StarTripleThree(Eigen::Vector3d conePoint1_, Eigen::Vector3d conePoint2_,Eigen::Vector3d conePoint3_,int col, int row,std::string imgPath)
{
    if(imgPath.size()==0) imgPath="/home/taiwei/cg/Assignment3/code_part/data/3333.jpg";

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
    this->conePoint3 = conePoint3_;
    this->conePoint1(0)+=(col-col_)/2;
    this->conePoint1(1)+=(row-row_)/2;
    this->conePoint2(0)+=(col-col_)/2;
    this->conePoint2(1)+=(row-row_)/2;
    this->conePoint3(0)+=(col-col_)/2;
    this->conePoint3(1)+=(row-row_)/2;

    this->row = row;
    this->col = col;
    angles[0]=M_PI/3*2;
    angles[1]=-M_PI/3*2;
    angles[2]=0;

    // cv::imshow("image",image);
    // cv::waitKey(0);
    // cv::destroyAllWindows();
    tilePiece = new TilePiece(image,row,col);
    mediumTilePiece = new TilePiece(image,row,col);
    complement1=new TilePiece(cv::Mat::zeros(tilePiece->row,tilePiece->col,CV_8UC3),row,col); //left
    complement2=new TilePiece(cv::Mat::zeros(tilePiece->row,tilePiece->col,CV_8UC3),row,col); //right
    complement3=new TilePiece(cv::Mat::zeros(tilePiece->row,tilePiece->col,CV_8UC3),row,col);
    //std::cout<<complement1->Coos[0].size()<<std::endl;
    //std::vector<Eigen::MatrixXd> tmp(CHANNEL_,Eigen::MatrixXd::Zero(row,col)),tmp_(CHANNEL_,Eigen::MatrixXd::Zero(row,col));
    std::vector<Eigen::MatrixXd> tmp,tmp_,tmp__,aux;

    //std::cout<<tmp[0].rows()<<" "<<tmp[0].cols()<<std::endl;
    int totalCol1 = 0, totalCol2 = 0, totalCol3 = 0;
    for(int i=0;i<2;i++)
    {
        aux=tilePiece->rotateCoo(conePoint1.head<2>(),angles[i]);
        tmp.insert(tmp.end(),aux.begin(),aux.end());
        //complement1->Coos.insert(complement1->Coos.end(),tmp.begin(),tmp.end());

        //std::cout<<complement1->Coos[i].rows()<<" "<<complement1->Coos[i].cols()<<std::endl;
        //HelperFunctions::concate(tmp,complement1->imgs,tilePiece->row,tilePiece->col,0);

        aux=tilePiece->rotateCoo(conePoint2.head<2>(),angles[i]);
        tmp_.insert(tmp_.end(),aux.begin(),aux.end());
        //complement2->Coos.insert(complement2->Coos.end(),tmp_.begin(),tmp_.end());
        //HelperFunctions::concate(tmp_,complement2->imgs,tilePiece->row,tilePiece->col,0);


        aux=tilePiece->rotateCoo(conePoint3.head<2>(),angles[i]);
        tmp__.insert(tmp__.end(),aux.begin(),aux.end());
    }
    // complement1->Coos = std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(3,totalCol1));
    // complement2->Coos = std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(3,totalCol2));
    // complement3->Coos = std::vector<Eigen::MatrixXd>(3,Eigen::MatrixXd::Zero(3,totalCol3));
    //for(int i=0;i<CHANNEL_;i+=1)
    //{
    complement1->Coos[0]=Eigen::MatrixXd::Zero(3,tmp[0].cols()+tmp[3].cols());
    complement1->Coos[1]=Eigen::MatrixXd::Zero(3,tmp[1].cols()+tmp[4].cols());
    complement1->Coos[2]=Eigen::MatrixXd::Zero(3,tmp[2].cols()+tmp[5].cols());
    complement2->Coos[0]=Eigen::MatrixXd::Zero(3,tmp_[0].cols()+tmp_[3].cols());
    complement2->Coos[1]=Eigen::MatrixXd::Zero(3,tmp_[1].cols()+tmp_[4].cols());
    complement2->Coos[2]=Eigen::MatrixXd::Zero(3,tmp_[2].cols()+tmp_[5].cols());
    complement3->Coos[0]=Eigen::MatrixXd::Zero(3,tmp__[0].cols()+tmp__[3].cols());
    complement3->Coos[1]=Eigen::MatrixXd::Zero(3,tmp__[1].cols()+tmp__[4].cols());
    complement3->Coos[2]=Eigen::MatrixXd::Zero(3,tmp__[2].cols()+tmp__[5].cols());
    complement1->Coos[0]<<tmp[0],tmp[3];
    complement1->Coos[1]<<tmp[1],tmp[4];
    complement1->Coos[2]<<tmp[2],tmp[5];
    complement2->Coos[0]<<tmp_[0],tmp_[3];
    complement2->Coos[1]<<tmp_[1],tmp_[4];
    complement2->Coos[2]<<tmp_[2],tmp_[5];
    complement3->Coos[0]<<tmp__[0],tmp__[3];
    complement3->Coos[1]<<tmp__[1],tmp__[4];
    complement3->Coos[2]<<tmp__[2],tmp__[5];
    //}
    // cv::Mat result;
    // HelperFunctions::cooToMat(complement1->Coos[0],tmp_[0],row,col);
    // HelperFunctions::cooToMat(complement1->Coos[2],tmp_[2],row,col);
    // HelperFunctions::cooToMat(complement1->Coos[1],tmp_[1],row,col);

    // std::vector<cv::Mat> results(3);
    // for(int i=0;i<3;i++)HelperFunctions::eigenToMat(tmp_[i],results[i]);
    // cv::merge(results,result);
    //
    // cv::imshow("result",result);
    // cv::waitKey();
    // cv::destroyAllWindows();

    std::cout<<"end"<<std::endl;
}
int StarTripleThree::cols(){return col;}
int StarTripleThree::rows(){return row;}



TilePiece* StarTripleThree::getMediumTilePiece(int LEVEL)
{

    int numInLayer;
    bool isOut = false;
    std::vector<Eigen::MatrixXd > tilesCoo =tilePiece->getCoos(), mTileCoo = mediumTilePiece->getCoos();
    Eigen::Vector3d cur,next;  //[2]>0: right, [2]<0: left
    Eigen::Vector2d disVec,disVec_ ;
    Eigen::Vector3d aux_;
    //double angles[]={M_PI/2,-M_PI/2,M_PI/2*3,0};
    std::vector<Eigen::Matrix3d> rotMats(3,Eigen::Matrix3d::Identity());
    rotMats[0].block(0,0,2,2) = HelperFunctions::getRotMatrix(angles[0],0);
    rotMats[1].block(0,0,2,2) = HelperFunctions::getRotMatrix(angles[1],0);
    bfs1.push(conePoint1);
    bfs1.push(conePoint2);
    bfs1.push(conePoint3);
    //Eigen::MatrixXd bfs1Coo(tilePiece->Coos[0]);
    std::vector<Eigen::MatrixXd> rotCoos(CHANNEL_),transCoos(CHANNEL_);

    //0-2: 1: red mouth;  3-5: 2: yellow mouth  6-8: 3: blue mouth
    aux_ = conePoint1- conePoint2;
    std::vector<Eigen::Vector3d> points2={rotMats[0]*aux_+conePoint2,rotMats[1]*aux_+conePoint2};

    aux_ = conePoint3- conePoint2;
    points2.push_back(rotMats[0]*aux_+conePoint2);
    points2.push_back(rotMats[1]*aux_+conePoint2);
    points2[0](2)=1;
    points2[1](2)=2;
    points2[2](2)=7;
    points2[3](2)=8;

    aux_ = conePoint2-conePoint3;
    std::vector<Eigen::Vector3d> points3={rotMats[0]*aux_+conePoint3,rotMats[1]*aux_+conePoint3};
    aux_ = conePoint1-conePoint3;
    points3.push_back(rotMats[0]*aux_+conePoint3);
    points3.push_back(rotMats[1]*aux_+conePoint3);
    points3[0](2)=4;
    points3[1](2)=5;
    points3[2](2)=1;
    points3[3](2)=2;

    aux_ = conePoint3 - conePoint1;
    std::vector<Eigen::Vector3d> points1={rotMats[0]*aux_+conePoint1,rotMats[1]*aux_+conePoint1};
    aux_ = conePoint2 - conePoint1;
    points1.push_back(rotMats[0]*aux_+conePoint1);
    points1.push_back(rotMats[1]*aux_+conePoint1);
    points1[0](2)=7;
    points1[1](2)=8;
    points1[2](2)=4;
    points1[3](2)=5;


    std::vector<std::vector<Eigen::MatrixXd>> complements1={complement1->Coos,complement1->rotateCoo(conePoint1.head<2>(),angles[0]),complement1->rotateCoo(conePoint1.head<2>(),angles[1])};
    std::vector<std::vector<Eigen::MatrixXd>> complements2={complement2->Coos,complement2->rotateCoo(conePoint2.head<2>(),angles[0]),complement2->rotateCoo(conePoint2.head<2>(),angles[1])};
    std::vector<std::vector<Eigen::MatrixXd>> complements3={complement3->Coos,complement3->rotateCoo(conePoint3.head<2>(),angles[0]),complement3->rotateCoo(conePoint3.head<2>(),angles[1])};

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
            if(cur[2]<3)
            {
                disVec = (cur-conePoint1).head<2>();
                //disVec_ = (cur-conePoint2).head<2>();
                angIdx = cur[2];
                for(int i=0;i<4;i++)
                {
                    idx = (angIdx+i)%4;
                    next<<points1[idx].head<2>()+disVec,points1[idx](2);
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
            else if(cur[2]>5)
            {
                //disVec_ = (cur-conePoint1).head<2>();
                disVec = (cur-conePoint3).head<2>();
                angIdx = cur[2]-6;
                for(int i=0;i<4;i++)
                {
                    idx=(angIdx+i)%4;
                    next<<points3[idx].head<2>()+disVec,points3[idx](2);
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
                HelperFunctions::transCoo(complements3[angIdx],disVec,transCoos);
                //std::cout<<"rot_"<<std::endl;
            }
            else
            {
                //disVec_ = (cur-conePoint1).head<2>();
                disVec = (cur-conePoint2).head<2>();
                angIdx = cur[2]-3;
                for(int i=0;i<4;i++)
                {
                    idx=(angIdx+i)%4;
                    next<<points2[idx].head<2>()+disVec,points2[idx](2);
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
        std::cout<<bfs1.size()<<" "<<bfs2.size()<<std::endl;
    }
    std::cout<<"bfs end"<<std::endl;
    return mediumTilePiece;
}
int StarTripleThree::visualize(int col,int row)
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

StarTripleThree::~StarTripleThree()
{
    std::cout<<"Done"<<std::endl;
}