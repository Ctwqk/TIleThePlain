//
// Created by taiwei on 11/30/24.
//

#include "HelperFunctions.h"

Eigen::MatrixXd HelperFunctions::coordinateMat; // Static Eigen matrix
bool HelperFunctions::isInited = false;
int HelperFunctions::CHANNEL_ = 3;
double HelperFunctions::little = 0.0000000001;
int HelperFunctions::init(int col,int row)
{
	coordinateMat = Eigen::MatrixXd(2,col*row);
	int idx=0;
	for(int i=1;i<=row;i++)
	{
		for(int j=1;j<=col;j++)
		{
			coordinateMat(0,idx) = (double)j;
			coordinateMat(1,idx) = (double)i;
			idx++;
		}
	}
	isInited = true;
	return 1;
}
Eigen::Matrix2d HelperFunctions::getRotMatrix(double angle,bool angleType)
{
	double degree;
	if(angleType)degree = angle * M_PI / 180;
	else degree = angle;
	double c=std::cos(degree);
	double s=std::sin(degree);
	Eigen::Matrix2d rotMat;
	rotMat<<c,-s,s,c;
	return rotMat;
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
int HelperFunctions::matToEigen(const cv::Mat &image, Eigen::MatrixXd &eigenImage ){
	if(image.empty()){
		std::cerr<<"expecting p0 to be not null"<<std::endl;
		return -1;
	}
	if(image.channels()!=1)
	{
		std::cerr<<"expacting 1 channel input"<<std::endl;
		return -1;
	}
	cv::Mat doubleImg;
	image.convertTo(doubleImg,CV_64F);
	//std::cout<<2<<std::endl;
	eigenImage=Eigen::Map<Eigen::MatrixXd>(doubleImg.ptr<double>(), doubleImg.rows,doubleImg.cols);

	//std::cout<<"2_"<<std::endl;
	return 0;
}
int HelperFunctions::eigenToMat(const Eigen::MatrixXd &eigenImage, cv::Mat &image)
{
	cv::Mat mat_(eigenImage.cols(),eigenImage.rows(),CV_64F,const_cast<double*>(eigenImage.data()));
	mat_.convertTo(image,CV_8U,255,0);
	return 0;
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
	//std::cout<<"get into"<<std::endl;
	for(int k=0;k<n;k++){
		//std::cout<<auxMat.size()<<" "<<eigenImages[k].size()<<std::endl;
		auxMat = 255.0*(eigenImages[k].array() - eigenImages[k].minCoeff()) / (eigenImages[k].maxCoeff()-eigenImages[k].minCoeff());
		//auxMat = eigenImages[k];

		for(int i=0;i<row;i++){
			rowPtr = channels[k].ptr<uchar>(i);
			for(int j=0;j<col;j++){
				rowPtr[j]=auxMat(i,j);
			}
		}
	}
	//std::cout<<"get out"<<std::endl;
    cv::merge(channels,image);
	return n;	
}

int HelperFunctions::rotateEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const Eigen::Vector2d &axis, Eigen::MatrixXd &eigenRotate, bool angleType){
	double degree ;
	if(angleType)degree = angle * M_PI / 180;
	else degree = angle;
	double c=std::cos(M_PI-degree);
	double s=std::sin(M_PI-degree);
	Eigen::Matrix2d rotMat;
	rotMat<<c,s,-s,c;
	double distance = std::hypot(axis[0],axis[1]);
	Eigen::MatrixXd aux1,aux2,aux3;
	//std::cout<<"a"<<std::endl;;
	//transEigen(eigenImage,-axis,aux1,1);
	aux1=eigenImage;
	Eigen::MatrixXd rotMapCoor = rotMat * HelperFunctions::coordinateMat;

	int row = eigenImage.rows();
	int col = eigenImage.cols();
	int x,y;
	//std::cout<<"a"<<std::endl;
	aux2 = Eigen::MatrixXd::Zero(row,col);
	for(int i=0;i<row*col;i++)
	{
		x=(int)(round(rotMapCoor(0,i))+col+axis[0]);
		y=(int)(round(rotMapCoor(1,i))+row+axis[1]);
		aux2((y%row),(x%col))=aux1((int)(i/col+row-axis[1])%row,(int)((i%col)+col-axis[0])%col);
	}
	//std::cout<<"a"<<std::endl;
	//transEigen(aux2,axis,aux3,1);
	eigenRotate = aux2;
	return 0;
}
int HelperFunctions::transEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const double &distance,Eigen::MatrixXd &eigenTrans, bool angleType, bool edgeOpe){
	//std::cout<<"in"<<std::endl;
	double degree;
	if(angleType) degree= angle * M_PI / 180;
	else degree = angle;
	int dx = round(distance * cos(degree));
	int dy = round(distance * sin(degree));
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
int HelperFunctions::transEigen(const Eigen::MatrixXd &eigenImage, Eigen::Vector2d disVec,Eigen::MatrixXd &eigenTrans, bool edgeOpe){
	int row = eigenImage.rows();
	int col = eigenImage.cols();
	int dx = round(disVec[0]);
	int dy = round(disVec[1]);
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
int HelperFunctions::reftEigen(const Eigen::MatrixXd &eigenImage, const double &angle, const Eigen::Vector2d &startPoint, Eigen::MatrixXd &eigenRefted, bool angleType , bool edgeOpe )
{
	double degree;
	if(angleType) degree= angle * M_PI / 180;
	else degree = angle;
	Eigen::MatrixXd aux1,aux2;
	rotateEigen(eigenImage,degree,startPoint,aux1,edgeOpe);
	aux1=aux1.rowwise().reverse();
	rotateEigen(aux1,degree,startPoint,aux2,edgeOpe);
	eigenRefted = aux2;
	return 0;
}
int HelperFunctions::matToCoo(const Eigen::MatrixXd &eigenImage, Eigen::MatrixXd &Coo)
{
	int row=eigenImage.rows();
	int col=eigenImage.cols();
	int idx=0;
	int n=row*col;
	//std::cout<<n<< std::endl;
	std::vector<double> x(n,-1),y(n,-1),value(n,-1);
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			value[idx]=eigenImage(i,j);
			if(value[idx]!=255)
			{
				x[idx]=j;
				y[idx]=i;
				idx++;
			}
		}
	}
	//std::cout<<idx<<std::endl;
	std::vector<double> data (x.begin(),x.begin()+idx);
	n=data.size();
	data.insert(data.end(),y.begin(),y.begin()+idx);
	data.insert(data.end(),value.begin(),value.begin()+idx);
	//std::cout<<1<<std::endl;
	Coo = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(data.data(),3,n);
	//std::cout<<"1_"<<std::endl;
	return 0;

}

int HelperFunctions::cooToMat(const Eigen::MatrixXd &Coo, Eigen::MatrixXd &mat,int row,int col,bool edgeOpe)
{
	mat = Eigen::MatrixXd::Zero(row,col);
	int x,y;
	for(int i=0;i<Coo.cols();i++)
	{
		y = Coo(1,i);
		x = Coo(0,i);
		if(x<0||y<0||x>=col||y>=row)
		{
			if(edgeOpe)
			{
				y=(y+row)%row;
				x=(x+col)%col;
			}
			else continue;
		}
		//std::cout<<x<<" "<<y<<" "<<col<<" "<<row<<std::endl;
		mat(y,x)=Coo(2,i);
	}
	//std::cout<<"in "<<row<<" "<<col<<" "<<mat.size()<<std::endl;
	//mat = mat_;
	return 0;
}
int HelperFunctions::normalize(const Eigen::MatrixXd &eigenImage, Eigen::MatrixXd &normalizedEigenImage)
{
	double minCoeff = eigenImage.minCoeff();
	double range = eigenImage.maxCoeff() - minCoeff;
	if(range==0) range=1;
	normalizedEigenImage = (eigenImage.array() - minCoeff) * (255.0 / range);
	return 0;
}
// int HelperFunctions::padding(const Eigen::MatrixXd &eigenImage, Eigen::MatrixXd &paddedEigenImage,int row,int col,int padVal)
// {
// 	if(row<eigenImage.rows()||col<eigenImage.cols())
// 	{
// 		std::cerr<<"expecting target size to be larger then tile size"<<std::endl;
// 		return -1;
// 	}
// 	Eigen::MatrixXd pad_ = Eigen::MatrixXd::Zero(row,col);
// 	//to be continue
// 	return 0;
// }
int HelperFunctions::concateCoo( std::vector<std::vector<Eigen::MatrixXd>> &Coos,std::vector<Eigen::MatrixXd>& results)
{
	Coos.push_back(results);
	int cols[3]{0,0,0};
	for(int i=0;i<Coos.size();i++)
	{
		cols[0]+=Coos[i][0].cols();
		cols[1]+=Coos[i][1].cols();
		cols[2]+=Coos[i][2].cols();
	}
	for(int i=0;i<results.size();i++)
	{
		results[i] = Eigen::MatrixXd::Zero(3,cols[i]);

	}
	cols[0]=0;
	cols[1]=0;
	cols[2]=0;
	for(int i=0;i<Coos.size();i++)
	{
		results[0].block(0,cols[0],3,Coos[i][0].cols()) = Coos[i][0];
		results[1].block(0,cols[1],3,Coos[i][0].cols()) = Coos[i][1];
		results[2].block(0,cols[2],3,Coos[i][0].cols()) = Coos[i][2];
		cols[0]+=Coos[i][0].cols();
		cols[1]+=Coos[i][1].cols();
		cols[2]+=Coos[i][2].cols();

	}
	Coos.pop_back();
	return 0;
}
int HelperFunctions::concate(const std::vector<Eigen::MatrixXd> &Coos, std::vector<Eigen::MatrixXd> &resultEigens, int row,int col, bool edgeOpe)
{
	int x,y,idx=0;
	for(const Eigen::MatrixXd &Coo : Coos){
		for(int i=0;i<Coo.cols();i++)
		{
			if(Coo(2,i)==0) continue;
			y = Coo(1,i);
			x = Coo(0,i);
			if(x<0||y<0||x>=col||y>=row)
			{
				if(edgeOpe)
				{
					y=(y+row)%row;
					x=(x+col)%col;
				}
				else continue;
			}
			//std::cout<<x<<" "<<y<<" "<<i<<" "<<resultEigens[idx].cols()<<" "<<resultEigens[idx].rows()<<" "<<col<<" "<<row<<std::endl;
			resultEigens[idx](y,x)=Coo(2,i);
		}
		idx++;
	}
	return 0;
}

int HelperFunctions::transCoo(std::vector<Eigen::MatrixXd>& imgCoos, Eigen::Vector2d disVec,  std::vector<Eigen::MatrixXd> &transCoos,bool edgeOpe)
{
	Eigen::Vector3d addOn;
	addOn << disVec, 0;
	//std::cout<<imgCoos.size()<<std::endl;
	for(int i=0;i<CHANNEL_;i++)
	{
		transCoos[i] = Eigen::MatrixXd::Zero(imgCoos[i].rows(),imgCoos[i].cols());
		//std::cout<<i<<std::endl;
		transCoos[i] = imgCoos[i].colwise() + addOn;
	}
	return 0;

}
int HelperFunctions::rotateCoo(std::vector<Eigen::MatrixXd> &imgCoos,  Eigen::Vector2d axis, double angle, std::vector<Eigen::MatrixXd> &rotateCoos, bool angleType , bool edgeOpe )
{
	std::vector<Eigen::MatrixXd> aux1,aux2;
	transCoo(imgCoos, -axis, rotateCoos, edgeOpe);
	Eigen::Matrix3d rotMat= Eigen::Matrix3d::Identity();
	rotMat.block(0,0,2,2)= HelperFunctions::getRotMatrix(angle, angleType);
	for(int i=0;i<CHANNEL_;i++)
	{
		rotateCoos[i]*=rotMat;
	}
	transCoo(imgCoos, axis, rotateCoos, edgeOpe);
	return 0;
}

GLuint HelperFunctions::matToTex(const cv::Mat&mat)
{
	GLuint textureID;
	glGenTextures(1, &textureID);
	glBindTexture(GL_TEXTURE_2D, textureID);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	GLenum format = GL_BGR;
	if (mat.channels() == 1)
		format = GL_RED;
	else if(mat.channels() == 3)
		format = GL_RGB;
	else if (mat.channels() == 4)
		format = GL_BGRA;
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, mat.cols, mat.rows, 0, format, GL_UNSIGNED_BYTE, mat.data);

	// Unbind the texture
	glBindTexture(GL_TEXTURE_2D, 0);

	return textureID;
}
int HelperFunctions::showImage(const cv::Mat&mat, int row, int col )
{
	if (!glfwInit())
	{
		std::cerr << "Failed to initialize glfw T^T";
		return -1;
	}
	GLFWwindow* window = glfwCreateWindow(1280, 720, "Simple Image Viewer", NULL, NULL);
    if (!window)
    {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130"); // Adjust GLSL version if needed


	GLuint textureID =matToTex(mat);

    while (!glfwWindowShouldClose(window))
    {
        glfwPollEvents();

        // Start ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Show the image in a simple ImGui window
        ImGui::Begin("Image Viewer");
        ImGui::Text("Displaying image:");
        ImGui::Image((ImTextureID)(intptr_t)textureID, ImVec2((float)col, (float)row));
        ImGui::End();

        // Rendering
        ImGui::Render();
        int displayW, displayH;
        glfwGetFramebufferSize(window, &displayW, &displayH);
        glViewport(0, 0, displayW, displayH);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    // Cleanup
    glDeleteTextures(1, &textureID);
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}