#include <iostream>
#include <vector>
#include <mpi.h>
#include <string>
#include "./TreeNode.h"
#include <math.h>
#include <algorithm>
#include <unistd.h>

void Direct_calc(const std::vector<body_t>& bodies, std::vector<double3_t>& forcedirect,int N)
{
	std::vector<double> force_Mag(bodies.size(),0);
	int count = 0;
	for ( int ibody = 0; ibody < N; ibody++)
	{
		forcedirect[count].r[0] = 0.0;
		forcedirect[count].r[1] = 0.0;
		forcedirect[count].r[2] = 0.0;
		
	for ( int n = 0; n < bodies.size(); n++)
	{
		if (ibody == n) continue;
		double Rx = bodies[ibody].r[0] - bodies[n].r[0];
		double Ry = bodies[ibody].r[1] - bodies[n].r[1];
		double Rz = bodies[ibody].r[2] - bodies[n].r[2];
		double R = sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
		forcedirect[count].r[0] += -G*bodies[ibody].m*bodies[n].m*Rx/(R*R*R);
		forcedirect[count].r[1] += -G*bodies[ibody].m*bodies[n].m*Ry/(R*R*R);
		forcedirect[count].r[2] += -G*bodies[ibody].m*bodies[n].m*Rz/(R*R*R);
	}
	force_Mag[count] = sqrt((forcedirect[count].r[0]*forcedirect[count].r[0])+ (forcedirect[count].r[1]*forcedirect[count].r[1]) + (forcedirect[count].r[2]*forcedirect[count].r[2]));
	count++;
	}
	
	
	//------------------------------------------------------------------
   	// Display FORCES STAT
   	//------------------------------------------------------------------
   	
	double DF_Mag_min,DF_Mag_max,DF_Mag_avg;
	DF_Mag_min= *std::min_element(std::begin(force_Mag),std::end(force_Mag));
	DF_Mag_max= *std::max_element(std::begin(force_Mag),std::end(force_Mag));
	DF_Mag_avg= std::accumulate(std::begin(force_Mag),std::end(force_Mag),0.0);
	DF_Mag_avg /=(double)(force_Mag.size());	
	
	std::cout<<"Maximum direct force : "<<DF_Mag_max<<std::endl;
	std::cout<<"Minimum direct force : "<<DF_Mag_min<<std::endl;
	std::cout<<"Average direct force : "<<DF_Mag_avg<<std::endl;
	
}


int main(int argc, char** argv)
{
    //----------------------------
    // MPI Setup
    //----------------------------

    int rank, nproc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
   if (argc < 4)
	{
		std::cerr << "You must pass :"<<std::endl<<"1.The number of bodies"<<std::endl<<"2. Theta "<<std::endl<<"3. Verify results (pass '1' to do verification)"<< std::endl;
		exit (1);
	}
	int N = atoi(argv[1]) ;
	double theta=strtod(argv[2],NULL);
	int Verify =atoi(argv[3]) ;
	srand(clock ()) ;
	
	
	
	//****************//
	//  SERIAL CASE   //
	//****************//
	if (nproc==1)
	{
		
		double t_start_BH = (double)clock()/(double)CLOCKS_PER_SEC;
		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
		// Point Creation (LOCAL POINTS)
		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
		
		std::vector<body_t> My_bodies(N);
		for (unsigned int n = 0; n < N; n++)
		{
			My_bodies[n].r[0] = (double)rand()/(double)RAND_MAX;
			My_bodies[n].r[1] = (double)rand()/(double)RAND_MAX;
			My_bodies[n].r[2] = (double)rand()/(double)RAND_MAX;
			My_bodies[n].m = (double)rand()/(double)RAND_MAX*1000;
		}
		double t_BH_1 = (double)clock()/(double)CLOCKS_PER_SEC;

		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
		// GET THE GLOBAL EXTENT
		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
		domain_t G_domain;
		double G_dimmin,G_dimmax;
		getPointExtent(My_bodies,G_domain,G_dimmin,G_dimmax,false);
		
		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
		// Build Tree (BASED ON THE GLOBAL EXTENT)
		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
		TreeNode root(G_dimmin-TOL,G_dimmax+TOL,G_dimmin-TOL,G_dimmax+TOL, G_dimmin -TOL, G_dimmax + TOL);
		
		for	 (unsigned int ibody = 0; ibody < N; ibody++)
		{
			root.addBody(My_bodies[ibody]);
		}

		double t_BH_2 = (double)clock()/(double)CLOCKS_PER_SEC;
		
/*
		int level = 0;
		int maxlevel = 0;	
		int nnodes = 0;
		int numberofbodies = 0;
		root.diagnostics (level , maxlevel ,numberofbodies, nnodes);
		std::cout<<"Tree Diagonistics: " << std::endl;
		std::cout<<"\tMax Level = " << maxlevel << std::endl;
		std::cout<<"\tN = " << numberofbodies << std::endl;
		std::cout<<"\tNode Count = " << nnodes << std::endl;
*/	
		root.computeCoM();
	
		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
		// Calculate the forces on my bodies (BH METHOD)
		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
		std::vector<double3_t> force_BH(My_bodies.size());
		std::vector<double> force_Mag(My_bodies.size(),0);
		std::vector<double> directforce_Mag(My_bodies.size(),0);
	
		
		for ( int i = 0; i < My_bodies.size(); i++)
		{
			root.computeForceOnBody(My_bodies[i],theta,force_BH[i]);
			force_Mag[i] = sqrt((force_BH[i].r[0]*force_BH[i].r[0])+ (force_BH[i].r[1]*force_BH[i].r[1]) + (force_BH[i].r[2]*force_BH[i].r[2]));
		}
		double t_stop_BH = (double)clock()/(double)CLOCKS_PER_SEC;
	
		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    	// Display FORCES STAT BH
    	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    	
		double F_Mag_min,F_Mag_max,F_Mag_avg;
		F_Mag_min= *std::min_element(std::begin(force_Mag),std::end(force_Mag));
		F_Mag_max= *std::max_element(std::begin(force_Mag),std::end(force_Mag));
		F_Mag_avg= std::accumulate(std::begin(force_Mag),std::end(force_Mag),0.0);
		F_Mag_avg /=(double)(force_Mag.size());
		
		std::cout<<"---------------------------------------"<<std::endl;
		std::cout<<"RESULTS (SERIAL BARNES-HUT METHOD)"<<std::endl;
		std::cout<<"Maximum force : "<<F_Mag_max<<std::endl;
		std::cout<<"Minimum force : "<<F_Mag_min<<std::endl;
		std::cout<<"Average force : "<<F_Mag_avg<<std::endl;	
		std::cout<<"Time (Overall time) : "<<t_stop_BH-t_start_BH<<" secs."<<std::endl;
		std::cout<<"Time (Point Generataion) : "<<t_BH_1-t_start_BH<<" secs."<<std::endl;
		std::cout<<"Time (Tree Building) : "<<t_BH_2-t_BH_1<<" secs."<<std::endl;
		std::cout<<"Time (Force computation) : "<<t_stop_BH-t_BH_2<<" secs."<<std::endl;
		std::cout<<"---------------------------------------"<<std::endl;
	
		if (Verify==1)
		{
			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    		// Calculate Forces Directly and Output FORCES STAT
    		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
			
			std::vector<double3_t> forcedirect(N);
			double t_start_direct = (double)clock()/(double)CLOCKS_PER_SEC;
			std::cout<<"---------------------------------------"<<std::endl;
			std::cout<<"RESULTS (SERIAL DIRECT METHOD)"<<std::endl;
			Direct_calc(My_bodies,forcedirect,N);
			double t_stop_direct = (double)clock()/(double)CLOCKS_PER_SEC;
			std::cout<<"Time (Overall) : "<<t_stop_direct-t_start_direct<<" secs."<<std::endl;
			std::cout<<"---------------------------------------"<<std::endl;
		}
		root.~TreeNode();
		MPI_Finalize();
		return 0;
	
	}
	
	
	//****************//
	// PARALLEL CASE  //
	//****************//
	
	double t_start_BH = (double)clock()/(double)CLOCKS_PER_SEC;
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Point Creation (LOCAL POINTS)
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	std::vector<body_t> bodies(N);
	for (unsigned int n = 0; n < N; n++)
	{
		bodies[n].r[0] = (double)rand()/(double)RAND_MAX;
		bodies[n].r[1] = (double)rand()/(double)RAND_MAX;
		bodies[n].r[2] = (double)rand()/(double)RAND_MAX;
		bodies[n].m = (double)rand()/(double)RAND_MAX*1000;
	} 
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// ORB and Get the points for my domain
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	//'send_count' count of pixels to be sent by this processor after image processing
	int send_count;
	//'recv_displacment' displacment of the data sent by this processor
	int recv_size;
	//'recv_counts' array hold the counts of data to be recieved on processor '0' from other processors (Gathering data after processing)
	int recv_counts[nproc];
	//'recv_disp' array hold the displacments in the receive buffer on processor '0' when gathering data from other processors
	int recv_disp[nproc];
	
    std::vector<std::vector<int> > points_in_orb_domains;
    //Partition the bodies on the processor to the domains
	ORB(nproc,bodies,points_in_orb_domains);
	//Vector of bodies to be sent to another processor
	std::vector<body_t> send_bodies(0);
	//Vector of bodies that belong to the domain assigned to this processor
	std::vector<body_t> My_bodies(0);
	
	for (int i=0;i<points_in_orb_domains.size();i++)
	{
		//Resize to 0 to add new elements to be sent to another processor
		send_bodies.resize(0);
		//Push the bodies in the domain 'i' to be sent to processor 'i'
		for (int k=0;k<points_in_orb_domains[i].size();k++)
		{
			send_bodies.push_back(bodies[points_in_orb_domains[i][k]]);
		}
		//Get the number of BYTES to be sent by this procesor to processor 'i'
		send_count=send_bodies.size()*sizeof(send_bodies[0]);
		//Gather in 'recv_counts' the counts to be sent by every processor to rank 'i' (store on rank 'i') 
    	MPI_Gather(&send_count, 1, MPI_INT, &recv_counts[0], 1, MPI_INT,i, MPI_COMM_WORLD);
		
		//Calculate the displacments to collect the bodies from all processors using a Gatherv call
		recv_disp[0]=0;
		for (int j=1;j<nproc;j++)
		{
			recv_disp[j]= recv_disp[j-1]+recv_counts[j-1];
		}
		
		//The total number of BYTES to be received
		recv_size=recv_disp[nproc-1]+recv_counts[nproc-1];
		
		if (rank==i)
		{	//Resize the recieve buffer on rank 'i'
			My_bodies.resize(recv_size/sizeof(send_bodies[0]));
		}
		//Gather the bodies that belong to the domain of processor 'i' in My_bodies 
		MPI_Gatherv(&send_bodies[0],send_count,MPI_BYTE,&My_bodies[0],recv_counts,recv_disp,MPI_BYTE,i,MPI_COMM_WORLD);
	}
	
	double t_BH_1 = (double)clock()/(double)CLOCKS_PER_SEC;
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// GET THE GLOBAL EXTENT
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	domain_t G_domain;
	double G_dimmin,G_dimmax;
	getPointExtent(My_bodies,G_domain,G_dimmin,G_dimmax,true);

	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Build Tree (BASED ON THE GLOBAL EXTENT)
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	TreeNode root(G_dimmin-TOL,G_dimmax+TOL,G_dimmin-TOL,G_dimmax+TOL, G_dimmin -TOL, G_dimmax + TOL);
		
	for	 (unsigned int ibody = 0; ibody < My_bodies.size(); ibody++)
	{
		root.addBody(My_bodies[ibody]);
	}
	root.computeCoM();
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// CALCULATE SOME DIAGNOSTICS
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	int level_b = 0;
	int maxlevel_b = 0;	
	int nnodes_b = 0;
	int numberofbodies_b = 0;
	body_t root_CoM_b;
	
	root.getCoM(root_CoM_b);
	root.diagnostics (level_b , maxlevel_b ,numberofbodies_b, nnodes_b);
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// GET THE LOCAL EXTENT OF ALL PROCESSORS
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	domain_t L_domain;
	double L_dimmin,L_dimmax;
	getPointExtent(My_bodies,L_domain,L_dimmin,L_dimmax,false);	
	
	//Gather local domains of all processors
	std::vector<domain_t> All_L_Domains(nproc);
	MPI_Allgather(&L_domain,sizeof(L_domain),MPI_BYTE,&All_L_Domains[0],sizeof(L_domain),MPI_BYTE,MPI_COMM_WORLD);


	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// GET THE LET FROM ALL PROCESSORS TO COMPLETE MY TREE
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	std::vector<std::vector<body_t>> LET_bodies(nproc,std::vector<body_t>(0));
	std::vector<std::vector<body_t>> bodies_to_add;
		
	for (int irank=0;irank<nproc;irank++)
	{
		if (irank!=rank)
		{
			root.LETBodies(All_L_Domains[irank], theta, LET_bodies[irank]);
		}	
	}
	MPI_Alltoall_vecvecT(LET_bodies,bodies_to_add);
	
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// ADD BODIES FROM ALL PROCESSORS TO COMPLETE MY TREE & UPDATE TREE COM
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	for (int irank=0;irank<nproc;irank++)
	{
		for (int i= 0; i<bodies_to_add[irank].size();i++)
		{
			root.addBody(bodies_to_add[irank][i]);
		}
	}
	root.computeCoM();
	
	
	double t_BH_2 = (double)clock()/(double)CLOCKS_PER_SEC;
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Calculate the forces on my bodies 
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	std::vector<double3_t> force_BH(My_bodies.size());
	std::vector<double> force_Mag(My_bodies.size(),0);
	
	for ( int i = 0; i < My_bodies.size(); i++)
	{
		root.computeForceOnBody(My_bodies[i],theta,force_BH[i]);
		force_Mag[i] = sqrt((force_BH[i].r[0]*force_BH[i].r[0])+ (force_BH[i].r[1]*force_BH[i].r[1]) + (force_BH[i].r[2]*force_BH[i].r[2]));
		
	}
	double t_stop_BH = (double)clock()/(double)CLOCKS_PER_SEC;
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// CALCULATE SOME DIAGNOSTICS
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−

	int level_a = 0;
	int maxlevel_a = 0;	
	int nnodes_a = 0;
	int numberofbodies_a = 0;
	body_t root_CoM_a;
	
	root.getCoM(root_CoM_a);
	root.diagnostics (level_a , maxlevel_a ,numberofbodies_a, nnodes_a);
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// OUTPUT DIAGNOSTICS
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	MPI_Barrier(MPI_COMM_WORLD);
	usleep(1000*rank);
	std::cout<<"Tree Diagonistics on rank: "<<rank << std::endl;
	std::cout<<"Before adding LET bodies received from other ranks"<< std::endl;
	std::cout<<"\tMax Level = " << maxlevel_b << std::endl;
	std::cout<<"\tN = " << numberofbodies_b << std::endl;
	std::cout<<"\tNode Count = " << nnodes_b << std::endl;
	std::cout<<"\tRoot CoM = (" << root_CoM_b.r[0]<<","<<root_CoM_b.r[1]<<","<<root_CoM_b.r[2]<<")"<< std::endl;
	
	std::cout<<"After adding LET bodies received from other ranks"<< std::endl;
	std::cout<<"\tMax Level = " << maxlevel_a << std::endl;
	std::cout<<"\tN = " << numberofbodies_a << std::endl;
	std::cout<<"\tNode Count = " << nnodes_a << std::endl;
	std::cout<<"\tRoot CoM = (" << root_CoM_a.r[0]<<","<<root_CoM_a.r[1]<<","<<root_CoM_a.r[2]<<")"<< std::endl;
	std::cout<<"---------------------------------------"<<std::endl;	
	

	
	
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// GET MAXIMUM, MINIMUM AND AVERAGE FORCE MAGNITUDE ON ALL PROCESSORS
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	double F_Mag_min,F_Mag_max,F_Mag_avg,G_F_Mag_min,G_F_Mag_max,G_F_Mag_avg;
	F_Mag_min= *std::min_element(std::begin(force_Mag),std::end(force_Mag));
	F_Mag_max= *std::max_element(std::begin(force_Mag),std::end(force_Mag));
	F_Mag_avg= std::accumulate(std::begin(force_Mag),std::end(force_Mag),0.0);
	F_Mag_avg /=(double)(force_Mag.size());
	
	MPI_Reduce(&F_Mag_min, &G_F_Mag_min, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
    MPI_Reduce(&F_Mag_max, &G_F_Mag_max, 1, MPI_DOUBLE, MPI_MAX, 0,MPI_COMM_WORLD);
	MPI_Reduce(&F_Mag_avg, &G_F_Mag_avg, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
	G_F_Mag_avg /=nproc;
	
	if (rank==0)
	{
		std::cout<<"---------------------------------------"<<std::endl;
		std::cout<<"RESULTS (PARALLEL BARNES-HUT METHOD)"<<std::endl;
		std::cout<<"Maximum force : "<<G_F_Mag_max<<std::endl;
		std::cout<<"Minimum force : "<<G_F_Mag_min<<std::endl;
		std::cout<<"Average force : "<<G_F_Mag_avg<<std::endl;
		std::cout<<"Time (Overall time) : "<<t_stop_BH-t_start_BH<<" secs."<<std::endl;
		std::cout<<"Time (Point Generation and ORB) : "<<t_BH_1-t_start_BH<<" secs."<<std::endl;
		std::cout<<"Time (Tree Building) : "<<t_BH_2-t_BH_1<<" secs."<<std::endl;
		std::cout<<"Time (Force computation) : "<<t_stop_BH-t_BH_2<<" secs."<<std::endl;
		std::cout<<"---------------------------------------"<<std::endl;
	}

		
	//****************//
	//Verification 	  //
	//****************//
	
	//I decided to test my algorthim by the most direct way
	//1. Collect all the bodies on rank 0
	//2. Calculate the force using one tree with all bodies
	//3. Calculate the force using the Direct Method
	
	if (Verify==1)
	{
		std::vector<body_t> AllBodies(N*nproc);
		MPI_Gather(&bodies[0],N*sizeof(bodies[0]),MPI_BYTE,&AllBodies[0],N*sizeof(bodies[0]),MPI_BYTE,0,MPI_COMM_WORLD);
		if (rank==0)
		{
			double T_start_BH = (double)clock()/(double)CLOCKS_PER_SEC;
			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
			// GET THE EXTENT
			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
		
			domain_t G_domain2;
			double G_dimmin2,G_dimmax2;
			getPointExtent(AllBodies,G_domain2,G_dimmin2,G_dimmax2,false);
			
			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
			// Build the Tree with all bodies
			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
		
			TreeNode root2(G_dimmin2-TOL,G_dimmax2+TOL,G_dimmin2-TOL,G_dimmax2+TOL, G_dimmin2 -TOL, G_dimmax2 + TOL);
			for	 (unsigned int ibody = 0; ibody < AllBodies.size(); ibody++)
			{
				root2.addBody(AllBodies[ibody]);
			}
				root2.computeCoM();
		
		
			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
			// Calculate the forces (BH METHOD)
			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
			std::vector<double3_t> force_BH(AllBodies.size());
			std::vector<double> force_Mag(AllBodies.size(),0);
			std::vector<double> directforce_Mag(AllBodies.size(),0);
		
			double T_BH_1 = (double)clock()/(double)CLOCKS_PER_SEC;
			for ( int i = 0; i < AllBodies.size(); i++)
			{
				root2.computeForceOnBody(AllBodies[i],theta,force_BH[i]);
				force_Mag[i] = sqrt((force_BH[i].r[0]*force_BH[i].r[0])+ (force_BH[i].r[1]*force_BH[i].r[1]) + (force_BH[i].r[2]*force_BH[i].r[2]));
			}
			double T_stop_BH = (double)clock()/(double)CLOCKS_PER_SEC;

			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
   		 	// Display FORCES STAT BH
    		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    		
			double F_Mag_min,F_Mag_max,F_Mag_avg;
			F_Mag_min= *std::min_element(std::begin(force_Mag),std::end(force_Mag));
			F_Mag_max= *std::max_element(std::begin(force_Mag),std::end(force_Mag));
			F_Mag_avg= std::accumulate(std::begin(force_Mag),std::end(force_Mag),0.0);
			F_Mag_avg /=(double)(force_Mag.size());
			std::cout<<"---------------------------------------"<<std::endl;	
			std::cout<<"RESULTS (SERIAL BARNES-HUT METHOD)"<<std::endl;
			std::cout<<"Maximum force : "<<F_Mag_max<<std::endl;
			std::cout<<"Minimum force : "<<F_Mag_min<<std::endl;
			std::cout<<"Average force : "<<F_Mag_avg<<std::endl;
			std::cout<<"Time (Overall time) : "<<T_stop_BH-T_start_BH<<" secs."<<std::endl;
			std::cout<<"Time (Force computation) : "<<T_BH_1-T_start_BH<<" secs."<<std::endl;
			std::cout<<"---------------------------------------"<<std::endl;

		
			//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    		// Calculate Forces Directly and Output FORCES STAT
    		//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
			std::vector<double3_t> forcedirect(N*nproc);
			std::cout<<"---------------------------------------"<<std::endl;
			std::cout<<"RESULTS (SERIAL DIRECT METHOD)"<<std::endl;
			double t_start_direct = (double)clock()/(double)CLOCKS_PER_SEC;
			Direct_calc(AllBodies,forcedirect,N*nproc);
			double t_stop_direct = (double)clock()/(double)CLOCKS_PER_SEC;
			std::cout<<"Time (Overall) : "<<t_stop_direct-t_start_direct<<" secs."<<std::endl;
			std::cout<<"---------------------------------------"<<std::endl;
		
			root2.~TreeNode();
		}
	}	
	root.~TreeNode();
	MPI_Finalize();
    return 0;
}

