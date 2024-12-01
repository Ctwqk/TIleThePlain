//
// Created by taiwei on 11/30/24.
//

#include "HelperFunctions.h"

Eigen::MatrixXd HelperFunctions::coordinateMat; // Static Eigen matrix
bool HelperFunctions::isInited = false;
double HelperFunctions::little = 0.0000000001;
int HelperFunctions::init(int col,int row)
{
	HelperFunctions::coordinateMat = Eigen::MatrixXd(2,col*row);
	int idx=0;
	for(int i=1;i<=row;i++)
	{
		for(int j=1;j<=col;j++)
		{
			HelperFunctions::coordinateMat(0,idx) = j;
			HelperFunctions::coordinateMat(1,idx) = i;
			idx++;
		}
	}
	HelperFunctions::isInited = true;
	return 1;
}

int HelperFunctions::matToEigen(const cv::Mat &image, std::vector<Eigen::MatrixXd> &eigenImages ){
  if(image.empty()){
    std::cerr<<"expecting p0 to be not null"<<std::endl;
    return -1;
  }
  std::vector<cv::Mat> channels;
  cv::split(image, channels);
  int n = channels.size();
  eigenImages = std::vector<Eigen::MatrixXd>(n,Eigen::MatrixXd::Zero(image.rows,image.cols));
  cv::split(image,channels);
	uchar* rowPtr ;
	for(int k=0;k<n;k++){
		for(int i=0;i<image.rows;i++){
			rowPtr = channels[k].ptr<uchar>(i);
			for(int j=0;j<image.cols;j++){
				eigenImages[k](i,j) = rowPtr[j];
			}
		}
	}
	return n;
}
int HelperFunctions::eigenToMat(const std::vector<Eigen::MatrixXd> &eigenImages, cv::Mat &image){
	int n=eigenImages.size();
	if(n==0) {
		std::cerr<<"expecting the eigen array to be non-emptr"<<std::endl;
		return -1;
	}
	int row = eigenImages[0].rows();
	int col = eigenImages[0].cols();
	std::vector<cv::Mat> channels(n,cv::Mat::zeros(row,col,CV_8U));
	uchar * rowPtr;
	Eigen::MatrixXd auxMat;
	for(int k=0;k<n;k++){
		auxMat = 255.0*(eigenImages[k].array() - eigenImages[k].minCoeff()) / (eigenImages[k].maxCoeff()-eigenImages[k].minCoeff());
		for(int i=0;i<row;i++){
			rowPtr = channels[k].ptr<uchar>(i);
			for(int j=0;j<col;j++){
				rowPtr[j]=auxMat(i,j);
			}
		}
	}
    cv::merge(channels,image);
	return n;	
}

int HelperFunctions::rotateEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const Eigen::Vector2d &axis, Eigen::MatrixXd &eigenRotate, bool angleType){
	double degree ;
	if(angleType)degree = angle * M_PI / 180;
	else degree = angle;
	double c=std::cos(degree);
	double s=std::sin(degree);
	Eigen::Matrix2d rotMat;
	rotMat<<c,s,-s,c;
	double distance = std::hypot(axis[0],axis[1]);
	Eigen::MatrixXd aux1,aux2,aux3;
	//std::cout<<"a"<<std::endl;;
	HelperFunctions::transEigen(eigenImage,std::atan2(axis[0],axis[1]),distance,aux1,1);

	Eigen::MatrixXd rotMapCoor = rotMat * HelperFunctions::coordinateMat;

	int row = eigenImage.rows();
	int col = eigenImage.cols();
	int x,y;
	//std::cout<<"a"<<std::endl;
	aux2 = Eigen::MatrixXd::Zero(row,col);
	for(int i=0;i<row*col;i++)
	{
		x=(int)(round(rotMapCoor(0,i))+col)%col;
		y=(int)(round(rotMapCoor(1,i))+row)%row;
		aux2(i/col,i%col)=aux1(y,x);
	}
	//std::cout<<"a"<<std::endl;
	HelperFunctions::transEigen(aux2,std::atan2(axis[0],axis[1]),-distance,aux3,0,1);
	eigenRotate = aux3;
	return 0;
}
int HelperFunctions::transEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const double &distance,Eigen::MatrixXd &eigenTrans, bool angleType, bool edgeOpe){
	//std::cout<<"in"<<std::endl;
	double degree;
	if(angleType) degree= angle * M_PI / 180;
	else degree = angle;
	int dx = round(distance * cos(degree));
	int dy = -round(distance * sin(degree));
	int row = eigenImage.rows();
	int col = eigenImage.cols();
	Eigen::MatrixXd aux1=Eigen::MatrixXd::Zero(row,col);
	double x1=0,x2=0,x3=row,x4=col,y1=0,y2=0,y3=row,y4=col;
	if(dx>0){
		x2+=dx;
		y4-=dx;
	}
	else{
		y2-=dx;
		x4+=dx;
	}
	if(dy>0){
		x1+=dy;
		y3-=dy;
	}
	else{
		y1-=dy;
		x3+=dy;
	}
	aux1.block(x1,x2,x3-x1,x4-x2)=eigenImage.block(y1,y2,y3-y1,y4-y2);
	if(edgeOpe)
	{
		if(x1==0)
		{
			aux1.block(x3,x2,row-x3,x4-x2)=eigenImage.block(0,y2,-dy,y4-y2);
			if(x2==0)
			{
				aux1.block(x1,x4,x3-x1,col-x4) = eigenImage.block(y1,0,y3-y1,-dx);
				aux1.block(x3,x4,row-x3,col-x4) = eigenImage.block(0,0,y1,y2);
			}
			else
			{
				aux1.block(x1,0,x3-x1,dx) = eigenImage.block(y1,y4,y3-y1,col-y4);
				aux1.block(x3,0,row-x3,x2) = eigenImage.block(0,y4,y1,col-y4);
			}
		}
		else
		{
			aux1.block(0,x2,dy,x4-x2)= eigenImage.block(y3,y2,row-y3,y4-y2);
			if(x2==0)
			{
				aux1.block(x1,x4,x3-x1,col-x4) = eigenImage.block(y1,0,y3-y1,-dx);
				aux1.block(0,x4,x1,col-x4) = eigenImage.block(y3,0,row-y3,y2);
			}
			else
			{
				aux1.block(x1,0,x3-x1,dx) = eigenImage.block(y1,y4,y3-y1,col-y4);
				aux1.block(0,0,x1,x2) = eigenImage.block(y3,y4,row-y3,col-y4);
			}
		}
	}
	eigenTrans = aux1;
	return 0;
}
