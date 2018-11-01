#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp> //for resizing
#include <vector>
#include <sstream>
#include <string>
#include <mpi.h>
#include <fstream>


#define MAX_SIZE 512
#define AFTER 10
#define BEFORE 11


//-----------------------------------------------------
// Usual Parallel Partitioning Code
//-----------------------------------------------------
void parallelRange(int globalstart, int globalstop, int irank, int nproc, int& localstart, int& localstop, int& localcount)
{
	int nvals = globalstop - globalstart + 1;
	int divisor = nvals/nproc;
	int remainder = nvals%nproc;
	int offset;
	if (irank < remainder) offset = irank;
	else offset = remainder;
    
	localstart = irank*divisor + globalstart + offset;
	localstop = localstart + divisor - 1;
	if (remainder > irank) localstop += 1;
	localcount = localstop - localstart + 1;
}
void colorPixelFromScalar(double f, cv::Vec3b& pixel)
{
	assert(f >= 0 - 1e-9 && f <= 1.0 + 1e-9);
	
	if (f < 0.03125) {pixel.val[2] = 59; pixel.val[1] = 76; pixel.val[0] = 192;}
	else if (f < 0.0625) {pixel.val[2] = 68; pixel.val[1] = 90; pixel.val[0] = 204;}
	else if (f < 0.09375) {pixel.val[2] = 77; pixel.val[1] = 104; pixel.val[0] = 215;}
	else if (f < 0.125) {pixel.val[2] = 87; pixel.val[1] = 117; pixel.val[0] = 225;}
	else if (f < 0.15625) {pixel.val[2] = 98; pixel.val[1] = 130; pixel.val[0] = 234;}
	else if (f < 0.1875) {pixel.val[2] = 108; pixel.val[1] = 142; pixel.val[0] = 241;}
	else if (f < 0.21875) {pixel.val[2] = 119; pixel.val[1] = 154; pixel.val[0] = 247;}
	else if (f < 0.25) {pixel.val[2] = 130; pixel.val[1] = 165; pixel.val[0] = 251;}
	else if (f < 0.28125) {pixel.val[2] = 141; pixel.val[1] = 176; pixel.val[0] = 254;}
	else if (f < 0.3125) {pixel.val[2] = 152; pixel.val[1] = 185; pixel.val[0] = 255;}
	else if (f < 0.34375) {pixel.val[2] = 163; pixel.val[1] = 194; pixel.val[0] = 255;}
	else if (f < 0.375) {pixel.val[2] = 174; pixel.val[1] = 201; pixel.val[0] = 253;}
	else if (f < 0.40625) {pixel.val[2] = 184; pixel.val[1] = 208; pixel.val[0] = 249;}
	else if (f < 0.4375) {pixel.val[2] = 194; pixel.val[1] = 213; pixel.val[0] = 244;}
	else if (f < 0.46875) {pixel.val[2] = 204; pixel.val[1] = 217; pixel.val[0] = 238;}
	else if (f < 0.5) {pixel.val[2] = 213; pixel.val[1] = 219; pixel.val[0] = 230;}
	else if (f < 0.53125) {pixel.val[2] = 221; pixel.val[1] = 221; pixel.val[0] = 221;}
	else if (f < 0.5625) {pixel.val[2] = 229; pixel.val[1] = 216; pixel.val[0] = 209;}
	else if (f < 0.59375) {pixel.val[2] = 236; pixel.val[1] = 211; pixel.val[0] = 197;}
	else if (f < 0.625) {pixel.val[2] = 241; pixel.val[1] = 204; pixel.val[0] = 185;}
	else if (f < 0.65625) {pixel.val[2] = 245; pixel.val[1] = 196; pixel.val[0] = 173;}
	else if (f < 0.6875) {pixel.val[2] = 247; pixel.val[1] = 187; pixel.val[0] = 160;}
	else if (f < 0.71875) {pixel.val[2] = 247; pixel.val[1] = 177; pixel.val[0] = 148;}
	else if (f < 0.75) {pixel.val[2] = 247; pixel.val[1] = 166; pixel.val[0] = 135;}
	else if (f < 0.78125) {pixel.val[2] = 244; pixel.val[1] = 154; pixel.val[0] = 123;}
	else if (f < 0.8125) {pixel.val[2] = 241; pixel.val[1] = 141; pixel.val[0] = 111;}
	else if (f < 0.84375) {pixel.val[2] = 236; pixel.val[1] = 127; pixel.val[0] = 99;}
	else if (f < 0.875) {pixel.val[2] = 229; pixel.val[1] = 112; pixel.val[0] = 88;}
	else if (f < 0.90625) {pixel.val[2] = 222; pixel.val[1] = 96; pixel.val[0] = 77;}
	else if (f < 0.9375) {pixel.val[2] = 213; pixel.val[1] = 80; pixel.val[0] = 66;}
	else if (f < 0.96875) {pixel.val[2] = 203; pixel.val[1] = 62; pixel.val[0] = 56;}
	else if (f < 1.0) {pixel.val[2] = 192; pixel.val[1] = 40; pixel.val[0] = 47;}
	else {pixel.val[2] = 180; pixel.val[1] = 4; pixel.val[0] = 38;}
}


void convertMatrixToImage(const std::vector<std::vector<double> >& matrix, cv::Mat& image)
{
    assert(matrix.size() == image.rows);
    for (int irow = 0; irow < image.rows; irow++)
    {
        assert(matrix[irow].size() == image.cols);
        for (int icol = 0; icol < image.cols; icol++)
        {
			//double value = matrix[irow][icol]/255;
			//colorPixelFromScalar(value,image.at<cv::Vec3b>(irow,icol));
			 
			image.at<uchar>(irow,icol)= matrix[irow][icol];
        }
    }
}

//-----------------------------------------------------------------------------------
//Put the rows of 'in' image starting from 'row_start' to 'row_stop' in 'out' vector 
//------------------------------------------------------------------------------------

void ToVec (const std::vector<std::vector<double> >& in, std::vector<double>& out, int row_start, int row_stop)
{
	//check that the first row is bigger than the last row and in the range of the image
	assert(row_start < row_stop && row_start >= 0 &&row_stop <= in.size());
	int index=0;
	for (int irow = row_start; irow < row_stop; irow++)
	{
		for (int icol = 0; icol < in[0].size(); icol++)
		{
			assert(index<out.size());
			out[index]=in[irow][icol];
			index++;		
		}
	}
}

void ToVecOfVec (const std::vector<double>& in, std::vector<std::vector<double> >& out)
{
	//check that the number of pixels in the image  equal the number of pixels in the vector
	assert(in.size()==(out.size()*out[0].size()));
	for (int irow = 0; irow < out.size(); irow++)
	{
		for (int icol = 0; icol < out[0].size(); icol++)
		{
			out[irow][icol]=in[icol+irow*out[0].size()];
		}
	}
}

void Generate_population (std::vector<std::vector<double> > &population,int Ghosts_before, int Ghosts_after)
{
	int ny=population.size();
	int nx=population[0].size();
	srand(clock());
    for (unsigned int iy = Ghosts_before; iy < ny-Ghosts_after; iy++)
    {
        for (unsigned int ix = 0; ix < nx; ix++)
        {
            //seed a 1/2 density of alive (just arbitrary really)
            int state = rand()%2;
            if (state == 0) population[iy][ix] = 255; //dead
            else population[iy][ix] = 0;   //alive
        }
    }
}


void New_Generation (std::vector<std::vector<double> > &population,int Ghosts_before, int Ghosts_after)
{
	int ny=population.size();
	int nx=population[0].size();
	
	std::vector<std::vector<double> > newpopulation = population;
	
	for (int iy = Ghosts_before; iy < ny-Ghosts_after; iy++)
		        {
		            for (int ix = 0; ix < nx; ix++)
		            {
		                int occupied_neighbours = 0;

		                for (int jy = iy - 1; jy <= iy + 1; jy++)
		                {
		                    for (int jx = ix - 1; jx <= ix + 1; jx++)
		                    {
		                        if (jx == ix && jy == iy) continue;
		                        
		                        int row = jy;
								//Not needed for the parallel case (I always include one before ghost and one after ghost)
								//Only needed for the serial case
		                        if (row == ny) row = 0;
		                        if (row == -1) row = ny-1;
                        
		                        int col = jx;
		                        if (col == nx) col = 0;
		                        if (col == -1) col = nx - 1;
                        
		                        if (population[row][col] == 0) occupied_neighbours++;
		                    }
		                }
            
		                if (population[iy][ix] == 0)   //alive
		                {
		                    if (occupied_neighbours <= 1 || occupied_neighbours >= 4) newpopulation[iy][ix] = 255; //dies
		                    if (occupied_neighbours == 2 || occupied_neighbours == 3) newpopulation[iy][ix] = 0; //same as population
		                }
		                else if (population[iy][ix] == 255) //dead
		                {
		                    if (occupied_neighbours == 3)
		                    {
		                        newpopulation[iy][ix] = 0; //reproduction
		                    }
		                }
		            }
		        }
		        population = newpopulation;
}

int main(int argc, char** argv)
{
	int rank, nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	//'recv_counts' array hold the counts of data to be recieved on processor '0' from other processors (Gathering data after processing)
	int recv_counts[nproc];
	//'recv_disp' array hold the displacments in the receive buffer on processor '0' when gathering data from other processors
	int recv_disp[nproc];
	//Vector to receive the processed data from all processors (using a Gatherv)
	std::vector<double> T_vec_all(0);
	
	
	int localstart,localstop,localcount;
	
	//Responsible for the global soln and process the TOP part of the grid
	if (rank==0)
	{
	
		assert(argc == 4);
		//-----------------------
		// Convert Command Line
		//-----------------------
				
		int Ny = atoi(argv[1]);
		int Nx = atoi(argv[2]);
    	int maxiter = atoi(argv[3]);
		
    	assert(Ny <= MAX_SIZE);
    	assert(Nx <= MAX_SIZE);
		
		
		if (nproc==1)
		{
			//****************//
    		//  SERIAL CASE  //
    		//***************//
		
			
		    //---------------------------------
		    // Generate the initial image
		    //---------------------------------
		    cv::Mat population(Ny, Nx, CV_8UC1);			
			std::vector<std::vector<double> > population_M((Ny), std::vector<double>(Nx,0.0));

			Generate_population(population_M,0,0);
			convertMatrixToImage(population_M,population);
		    cv::namedWindow("Population", cv::WINDOW_AUTOSIZE );
		    cv::Mat image_for_viewing(MAX_SIZE,MAX_SIZE,CV_8UC1);
			
			double t_start = MPI_Wtime();
		    for (int iter = 0; iter < maxiter; iter++)
		    {
		        //something new here - we will resize our image up to MAX_SIZE x MAX_SIZE so its not really tiny on the screen
		        cv::resize(population,image_for_viewing,image_for_viewing.size(),cv::INTER_LINEAR); //resize image to MAX_SIZE x MAX_SIZE
		        cv::imshow("Population", image_for_viewing);
		        cvWaitKey(10);	//wait 10 seconds before closing image (or a keypress to close)
				New_Generation(population_M,0,0);
				convertMatrixToImage(population_M,population);
				
				//Save to a file every 200 iteration
				if(iter%200==0)
				{
					//Set the name of the output file for full solutions
					std::ostringstream name;
					name<< "(Full) At iteration: '"<<iter<<"' on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
					//Write the image to the file
					cv::imwrite(name.str(),population);
					
				}
		    }
			double t_stop = MPI_Wtime();
			
			//Set the name of the output file for full solutions
			std::ostringstream name;
			name<< "(Full) Final state on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
			//Write the image to the file
			cv::imwrite(name.str(),population);
			
			std::cout<<"Serial execution time for '" <<Nx*Ny<<"' points is '"<<t_stop-t_start<<"' secs."<<std::endl;
			MPI_Finalize();
			return 0;
		}
		
		
		//****************//
    	// PARALLEL CASE //
    	//***************//
		
		

		//----------------------------------------
    	// Global solution Space Creation
    	//----------------------------------------
		
		//Resize to receive data
		T_vec_all.resize(Nx*Ny);
		//To store processed data from all processors in matrix format
		std::vector<std::vector<double> > T_G(Ny, std::vector<double>(Nx,0.0));
		
		
		//Calculate 'recv_counts' & 'recv_disp' for Variable Gathering
		for(int irank=0;irank<nproc;irank++)
		{
			parallelRange(0, Ny-1, irank, nproc, localstart, localstop, localcount);
			recv_counts[irank]=Nx*localcount;
			recv_disp[irank]=Nx*localstart;
		}
		
		//----------------------------------------
    	// Local part soln space Solution Space Creation
    	//----------------------------------------
		
		int ny;
		//Get range of rows to work on
		parallelRange(0, Ny-1, rank, nproc, localstart, localstop, localcount);
		//Vector of processed points to be sent to processor 0
		std::vector<double> T_vec(localcount*Nx,0);
    	
    	//Two Ghost rows (After rank 0 part) (Before rank 0 from the last processor)
		ny=localcount+2;
		
		//Set Local Matrix to process my part
    	std::vector<std::vector<double> > T((ny), std::vector<double>(Nx,0.0));
    	//std::vector<std::vector<double> > Tnew = T;
		
		//Generate the population
		Generate_population(T,1,1);
		
		
		
		//----------------------------------------
    	// Setup Image for Display 
    	//----------------------------------------
		
    	cv::Mat image(Ny,Nx,CV_8UC1);
    	cv::Mat image_to_view(MAX_SIZE, MAX_SIZE, CV_8UC1);
    	cv::namedWindow("Population", cv::WINDOW_AUTOSIZE );
    	cv::imshow("Population", image_to_view);
		cvWaitKey(10);
		
		//----------------------------------------
    	//To Save the Local part every 200 iterations to a file  (Including Ghost rows)
    	//----------------------------------------

		cv::Mat image_P(ny,Nx,CV_8UC1);
    	cv::Mat image_to_view_P(MAX_SIZE/nproc, MAX_SIZE, CV_8UC1);
    	convertMatrixToImage(T,image_P);
		
		
		MPI_Request request;
		
		MPI_Barrier(MPI_COMM_WORLD); 
		double t_start = MPI_Wtime();
		
		for (int iter = 0; iter < maxiter; iter++)
		{
			//Send the first row I processed to the last rank
			MPI_Isend(&T[1][0],Nx, MPI_DOUBLE, nproc-1, AFTER, MPI_COMM_WORLD,&request);
			//Send the last row I processed to the next rank
			MPI_Isend(&T[ny-2][0],Nx, MPI_DOUBLE, rank+1, BEFORE, MPI_COMM_WORLD,&request);
			//Receive the (After) ghost row from the next rank
			MPI_Recv(&T[ny-1][0],Nx,MPI_DOUBLE,rank+1,AFTER,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			//Receive the (Before) ghost row from the last rank
			MPI_Recv(&T[0][0],Nx,MPI_DOUBLE,nproc-1,BEFORE,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			//Move to new generation
			New_Generation(T,1,1);
			
			//Output every 10 iterations
			//if (iter%10 == 0)
        	{
				//Convert my part from matrix to vector to do gathering
				ToVec (T, T_vec, 0, localcount);
				//Gather the processed data from all processor and store in 'T_vec_all' on  processor '0'
				MPI_Gatherv(&T_vec[0],localcount*Nx,MPI_DOUBLE,&T_vec_all[0],recv_counts,recv_disp,MPI_DOUBLE,0,MPI_COMM_WORLD);
				//Convert the gathered vector to a matrix (Grid)
				ToVecOfVec(T_vec_all,T_G);
				//Convert the matrix to image and display it
				convertMatrixToImage(T_G,image);
            	cv::resize(image,image_to_view,image_to_view.size(), cv::INTER_LINEAR);
            	cv::imshow("Population", image_to_view);
            	cvWaitKey(10);
				
        	}
			
			//Save to a file every 200 iteration
			if(iter%200==0)
			{
				convertMatrixToImage(T,image_P);
				//Set the name of the output file
				std::ostringstream name;
				name<< "At iteration: '"<<iter<<"' on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
				//Write the image to the file
				cv::imwrite(name.str(),image_P);
				//Set the name of the output file for full solutions
				std::ostringstream name2;
				name2<< "(Full) At iteration: '"<<iter<<"' on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
				//Write the image to the file
				cv::imwrite(name2.str(),image);
		
			}
		}
		MPI_Barrier(MPI_COMM_WORLD); 
		double t_stop = MPI_Wtime();
		
		convertMatrixToImage(T,image_P);
		//Set the name of the output file
		std::ostringstream name;
		name<< "Final state on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
		//Write the image to the file
		cv::imwrite(name.str(),image_P);
		
		//Set the name of the output file for full solutions
		std::ostringstream name2;
		name2<< "(Full) Final state on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
		//Write the image to the file
		cv::imwrite(name2.str(),image);
		
		std::cout<<"Parallel execution time for '" <<Nx*Ny<<"' points using '"<<nproc<< "' processors is '"<<t_stop-t_start<<"' secs."<<std::endl;
		
    
		
	}
	//Responsible for the BOTTOM Part of the grid
	else if(rank==nproc-1)
	{
		//-----------------------
		// Convert Command Line
		//-----------------------
				
		int Ny = atoi(argv[1]);
		int Nx = atoi(argv[2]);
    	int maxiter = atoi(argv[3]);
		
		//----------------------------------------
    	// Local part soln space Solution Space Creation
    	//----------------------------------------
		
		int ny;
		//Get range of rows to work on
		parallelRange(0, Ny-1, rank, nproc, localstart, localstop, localcount);
		//Vector of processed points to be sent to processor 0
		std::vector<double> T_vec(localcount*Nx,0);
    	
    	//Two Ghost row (Before part for the last processors) (After from rank 0) 
		ny=localcount+2;
		
		//Set Local Matrix to process my part
    	std::vector<std::vector<double> > T((ny), std::vector<double>(Nx,0.0));
    	//std::vector<std::vector<double> > Tnew = T;
		
		//Generate the population
		Generate_population(T,1,1);
		
		//----------------------------------------
    	//To Save the Local part every 200 iterations to a file  (Including Ghost rows)
    	//----------------------------------------

		cv::Mat image_P(ny,Nx,CV_8UC1);
    	cv::Mat image_to_view_P(MAX_SIZE/nproc, MAX_SIZE, CV_8UC1);
    	convertMatrixToImage(T,image_P);
		
		
		MPI_Request request;
		
		MPI_Barrier(MPI_COMM_WORLD);
		for (int iter = 0; iter < maxiter; iter++)
		{
			//Send the first row I processed to the previous rank
			MPI_Isend(&T[1][0],Nx, MPI_DOUBLE, rank-1, AFTER, MPI_COMM_WORLD,&request);
			//Send the last row I processed to rank 0
			MPI_Isend(&T[ny-2][0],Nx, MPI_DOUBLE, 0, BEFORE, MPI_COMM_WORLD,&request);
			//Receive the (Before) ghost row from the previous rank
			MPI_Recv(&T[0][0],Nx,MPI_DOUBLE,rank-1,BEFORE,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			//Receive the (After) ghost row from rank 0
			MPI_Recv(&T[ny-1][0],Nx,MPI_DOUBLE,0,AFTER,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			//Move to new generation
			New_Generation(T,1,1);
			
			//Output every 10 iterations
			//if (iter%10 == 0)
        	{
				//Convert my part from matrix to vector to do gathering
				ToVec (T, T_vec, 1, localcount+1);
				//Gather the processed data from all processor and store in 'T_vec_all' on  processor '0'
				MPI_Gatherv(&T_vec[0],localcount*Nx,MPI_DOUBLE,&T_vec_all[0],recv_counts,recv_disp,MPI_DOUBLE,0,MPI_COMM_WORLD);
        	}
			
			//Save to a file every 100 iteration
			if(iter%200==0)
			{
				convertMatrixToImage(T,image_P);
				//Set the name of the output file
				std::ostringstream name;
				name<< "At iteration: '"<<iter<<"' on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
				//Write the image to the file
				cv::imwrite(name.str(),image_P);
				
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		convertMatrixToImage(T,image_P);
		//Set the name of the output file
		std::ostringstream name;
		name<< "Final state on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
		//Write the image to the file
		cv::imwrite(name.str(),image_P);
		

	}
	//Middle parts
	else
	{
		//-----------------------
		// Convert Command Line
		//-----------------------
				
		int Ny = atoi(argv[1]);
		int Nx = atoi(argv[2]);
    	int maxiter = atoi(argv[3]);
		
		//----------------------------------------
    	// Local part soln space Solution Space Creation
    	//----------------------------------------
		
		int ny;
		//Get range of rows to work on
		parallelRange(0, Ny-1, rank, nproc, localstart, localstop, localcount);
		//Vector of processed points to be sent to processor 0
		std::vector<double> T_vec(localcount*Nx,0);
    	
    	//One Ghost row (After rank 0 part)
		ny=localcount+2;
		
		//Set Local Matrix to process my part
    	std::vector<std::vector<double> > T((ny), std::vector<double>(Nx,0.0));
    	//std::vector<std::vector<double> > Tnew = T;
		
		//Generate the population
		Generate_population(T,1,1);
		
		
		//----------------------------------------
    	//To Save the Local part every 200 iterations to a file  (Including Ghost rows)
    	//----------------------------------------

		cv::Mat image_P(ny,Nx,CV_8UC1);
    	cv::Mat image_to_view_P(MAX_SIZE/nproc, MAX_SIZE, CV_8UC1);
    	convertMatrixToImage(T,image_P);
		

		
		MPI_Request request;
		
		MPI_Barrier(MPI_COMM_WORLD);
		for (int iter = 0; iter < maxiter; iter++)
		{
			//Send the last row I processed to the next rank
			MPI_Isend(&T[ny-2][0],Nx, MPI_DOUBLE, rank+1, BEFORE, MPI_COMM_WORLD,&request);
			//Send the first row I processed to the previous rank
			MPI_Isend(&T[1][0],Nx, MPI_DOUBLE, rank-1, AFTER, MPI_COMM_WORLD,&request);
			//Receive the (After) ghost row from the next rank
			MPI_Recv(&T[ny-1][0],Nx,MPI_DOUBLE,rank+1,AFTER,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			//Receive the (Before) ghost row from the previous rank
			MPI_Recv(&T[0][0],Nx,MPI_DOUBLE,rank-1,BEFORE,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        	
			//Move to new generation
			New_Generation(T,1,1);
			
			//Output every 10 iterations
			//if (iter%10 == 0)
        	{
				//Convert my part from matrix to vector to do gathering
				ToVec (T, T_vec, 1, localcount+1);
				//Gather the processed data from all processor and store in 'T_vec_all' on  processor '0'
				MPI_Gatherv(&T_vec[0],localcount*Nx,MPI_DOUBLE,&T_vec_all[0],recv_counts,recv_disp,MPI_DOUBLE,0,MPI_COMM_WORLD);
        	}
			
			//Save to a file every 200 iteration
			if(iter%200==0)
			{
				convertMatrixToImage(T,image_P);
				//Set the name of the output file
				std::ostringstream name;
				name<< "At iteration: '"<<iter<<"' on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
				//Write the image to the file
				cv::imwrite(name.str(),image_P);
			
			}
			
		}
		MPI_Barrier(MPI_COMM_WORLD);
		
		convertMatrixToImage(T,image_P);
		//Set the name of the output file
		std::ostringstream name;
		name<< "Final state on processor '"<<rank<<"' of '"<<nproc<<"' processors.jpg";
		//Write the image to the file
		cv::imwrite(name.str(),image_P);
		

	}
	
	
	
	
	MPI_Finalize();
}