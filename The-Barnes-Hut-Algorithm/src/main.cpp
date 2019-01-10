//-------------------------------------------------
// Parallel Implementation of the Barnes-Hut algorithm using OpenMPI Library
// Course: 'ECE4530-Parallel Processing' at the University of Manitoba.
// Implemented by: Micheal Sobhy
// (The initial code structure and some functions were provided by the course instructor)
// Email: sobhymich@gmail.com
//----------------------------------------------------

#include <iostream>
#include <vector>
#include <mpi.h>
#include <string>
#include "./TreeNode.h"
#include <math.h>
#include <algorithm>
#include <unistd.h>

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Display results (maximum, minimum, and average force magnitude) and computational time
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
void DisplayForceStat(std::vector<double>& force_Mag, double t_BodiesGeneration,
					 double t_TreeBuilding, double t_ForceComp, bool global)
{	
	//−−−−−−−−−−−−−
    // MPI setup 
    //−−−−−−−−−−−−−
	int rank, nproc;	// Rank of the running processor, Number of running processors
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);	//Get the number of processors available
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);	//Identify the ID of the processor
	
	double F_Mag_min,F_Mag_max,F_Mag_avg,G_F_Mag_min,G_F_Mag_max,G_F_Mag_avg;
	F_Mag_min= *std::min_element(std::begin(force_Mag),std::end(force_Mag));
	F_Mag_max= *std::max_element(std::begin(force_Mag),std::end(force_Mag));
	F_Mag_avg= std::accumulate(std::begin(force_Mag),std::end(force_Mag),0.0);
	F_Mag_avg /=(double)(force_Mag.size());
	
	//−−−−−−−−−−−−−
    // Obtain Global results
    //−−−−−−−−−−−−−
	if (nproc !=1 && global)
	{	
		MPI_Reduce(&F_Mag_min, &G_F_Mag_min, 1, MPI_DOUBLE, MPI_MIN, 0,MPI_COMM_WORLD);
    	MPI_Reduce(&F_Mag_max, &G_F_Mag_max, 1, MPI_DOUBLE, MPI_MAX, 0,MPI_COMM_WORLD);
		MPI_Reduce(&F_Mag_avg, &G_F_Mag_avg, 1, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);
		G_F_Mag_avg /=nproc;
	}
	else
	{	
		G_F_Mag_min=F_Mag_min;
		G_F_Mag_max=F_Mag_max;
		G_F_Mag_avg=F_Mag_avg;
		
	}
	if (rank==0)
	{	
		std::cout<<"Maximum force : "<<G_F_Mag_max<<std::endl;
		std::cout<<"Minimum force : "<<G_F_Mag_min<<std::endl;
		std::cout<<"Average force : "<<G_F_Mag_avg<<std::endl<<std::endl;	
		
		if(t_BodiesGeneration!=0)
			std::cout<<"Time (Point Generataion) : "<<t_BodiesGeneration<<" secs."<<std::endl;
		if(t_TreeBuilding!=0)
			std::cout<<"Time (Tree Building) : "<<t_TreeBuilding<<" secs."<<std::endl;
		if(t_ForceComp!=0)
			std::cout<<"Time (Force computation) : "<<t_ForceComp<<" secs."<<std::endl;
		std::cout<<"Time (Overall time) : "<<t_BodiesGeneration+t_TreeBuilding+t_ForceComp<<" secs."<<std::endl;
	}
	
}

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Direct solution for the N-body Problem (to verify the results from the Barnes-Hut algorithm)
// Calculate the force on each body by summing the forces due to all other bodies
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
void Direct_calc(const std::vector<body_t>& bodies)
{
	std::vector<double3_t> forcedirect(bodies.size());
	std::vector<double> force_Mag(bodies.size(),0);
	int count = 0;
		
	//Start Timing
	double t_2 = (double)clock()/(double)CLOCKS_PER_SEC;
	for ( int ibody = 0; ibody < bodies.size(); ibody++)
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
	//Stop Timing 
	double t_3 = (double)clock()/(double)CLOCKS_PER_SEC;
	
	//−−−−−−−−−−−−−−−−−−−−−−
   	// Display results
   	//−−−−−−−−−−−−−−−−−−−−−−
	std::cout<<"---------------------------------------"<<std::endl;
	std::cout<<"RESULTS (DIRECT METHOD)"<<std::endl;
	DisplayForceStat(force_Mag,0,0,t_3-t_2,false);

}

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Generate N bodies at random positions over a unit cube, having random weights
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
void Generate_Bodies(std::vector<body_t>& bodies)
{
	srand(clock ()) ;	// Seed the random number generator
	for (unsigned int n = 0; n < bodies.size(); n++)
	{
		bodies[n].r[0] = (double)rand()/(double)RAND_MAX;
		bodies[n].r[1] = (double)rand()/(double)RAND_MAX;
		bodies[n].r[2] = (double)rand()/(double)RAND_MAX;
		bodies[n].m = (double)rand()/(double)RAND_MAX*1000;
	}
	
}

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Add the bodies generated to the tree
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
TreeNode* Build_Tree(std::vector<body_t>& My_bodies, double dimmin, double dimmax)
{
	
	// Create the root tree node with the global extent
	TreeNode* root= new TreeNode(dimmin-TOL,dimmax+TOL,dimmin-TOL,dimmax+TOL, dimmin -TOL, dimmax + TOL);
		
	//Add all the bodies to the tree
	for	 (unsigned int ibody = 0; ibody < My_bodies.size(); ibody++)
	{
		root->addBody(My_bodies[ibody]);
	}
	
	//Calculate the Center of Mass of all the tree nodes
	root->computeCoM();
	
	return root;
}

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Serial implementation of Barnes-Hut algorithm to calculate
// the gravitational forces between the generated bodies
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−	
void Serial_BH(std::vector<body_t>& My_bodies,double theta,int Verify,double t_BodiesGeneration)
{
	
	//−−−−−−−−−−−−
	// Build Tree 
	//−−−−−−−−−−−−
	
	double t_1 = (double)clock()/(double)CLOCKS_PER_SEC;
	// GET THE GLOBAL EXTENT (Domain that contain all the bodies)	
	domain_t G_domain;
	double G_dimmin,G_dimmax;
	getPointExtent(My_bodies,G_domain,G_dimmin,G_dimmax,false);
	
	TreeNode* root=Build_Tree(My_bodies,G_dimmin,G_dimmax);
	
	//Record time for building the tree
	double t_2 = (double)clock()/(double)CLOCKS_PER_SEC;
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Calculate the forces on all the bodies
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	std::vector<double3_t> force_BH(My_bodies.size());
	std::vector<double> force_Mag(My_bodies.size(),0);
	std::vector<double> directforce_Mag(My_bodies.size(),0);
	
	// Calculate the force of the N bodies using Barnes-Hut algorithm	
	for ( int i = 0; i < My_bodies.size(); i++)
	{
		root->computeForceOnBody(My_bodies[i],theta,force_BH[i]);
		force_Mag[i] = sqrt((force_BH[i].r[0]*force_BH[i].r[0])+ (force_BH[i].r[1]*force_BH[i].r[1]) + (force_BH[i].r[2]*force_BH[i].r[2]));
	}
	double t_3 = (double)clock()/(double)CLOCKS_PER_SEC;
	
	std::cout<<"---------------------------------------"<<std::endl;
	std::cout<<"RESULTS (SERIAL BARNES-HUT METHOD)"<<std::endl;
	DisplayForceStat(force_Mag,t_BodiesGeneration,t_2-t_1,t_3-t_2,false);
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
    // Verification using direct method
    //−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	if (Verify==1)
	{
		
		Direct_calc(My_bodies);

	}
	root->prune();
	root->~TreeNode();
}


//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Find the bodies needed by any other processor to complete its local tree, 
// then communicate the bodies to complete the local trees
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−	
void CompleteLocalTrees(TreeNode* root,std::vector<body_t>& My_bodies, double theta)
{	
	//−−−−−−−−−−−−−
    // MPI setup 
    //−−−−−−−−−−−−−
	int rank, nproc;	// Rank of the running processor, Number of running processors
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);	//Get the number of processors available
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);	//Identify the ID of the processor
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Get the local extent domains of all processors
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	domain_t L_domain;
	double L_dimmin,L_dimmax;
	getPointExtent(My_bodies,L_domain,L_dimmin,L_dimmax,false);	
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Gather local extent domains of all processors
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	std::vector<domain_t> All_L_Domains(nproc);
	MPI_Allgather(&L_domain,sizeof(L_domain),MPI_BYTE,&All_L_Domains[0],sizeof(L_domain),MPI_BYTE,MPI_COMM_WORLD);


	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Get the locally essential bodies to complete my local tree
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−	
	std::vector<std::vector<body_t>> LET_bodies(nproc,std::vector<body_t>(0));
	std::vector<std::vector<body_t>> bodies_to_add;
		
	for (int irank=0;irank<nproc;irank++)
	{
		if (irank!=rank)
		{
			root->LETBodies(All_L_Domains[irank], theta, LET_bodies[irank]);
		}	
	}
	// Communicate the bodies to other processors and get the bodies need to complete my local tree
	MPI_Alltoall_vecvecT(LET_bodies,bodies_to_add);
	
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Add bodies from all processors to complete my local tree then update COM
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	
	for (int irank=0;irank<nproc;irank++)
	{
		for (int i= 0; i<bodies_to_add[irank].size();i++)
		{
			root->addBody(bodies_to_add[irank][i]);
		}
	}
	root->computeCoM();
}

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Communicate the bodies such that each processor has all the bodies 
// that belong to the domain it's assigned
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
void CommunicateDomainsBodies( std::vector<body_t>& bodies,
							   std::vector<std::vector<int> >& points_in_orb_domains,
							   std::vector<body_t>& My_bodies)
{
	//−−−−−−−−−−−−−
    // MPI setup 
    //−−−−−−−−−−−−−
	int rank, nproc;	// Rank of the running processor, Number of running processors
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);	//Get the number of processors available
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);	//Identify the ID of the processor
	
	
	int send_count; // Number of BYTES to be sent by this processor
	int recv_size; // Total number of BYTES to be received
	int recv_counts[nproc];// Number of BYTES to be recieved from each processor
	int recv_disp[nproc];	// Displacments in the receive buffer on the current processor
	std::vector<body_t> send_bodies(0);	//Vector of bodies to be sent to another processor

	
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
	
	
}

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Parallel implementation of Barnes-Hut algorithm to calculate
// the gravitational forces between the generated bodies
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−	
void Parallel_BH(std::vector<body_t>& bodies,int N,double theta,int Verify,double t_BodiesGeneration)
{
    //−−−−−−−−−−−−−
    // MPI setup 
    //−−−−−−−−−−−−−
	int rank, nproc;	// Rank of the running processor, Number of running processors
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);	//Get the number of processors available
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);	//Identify the ID of the processor
	

	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Partition the bodies on this processor to orthogonal domains equal to the number of processors,
	// and communicate the bodies such that each processor is responsible for a domain
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
  	std::vector<std::vector<int> > bodies_in_orb_domains; // Reference to the bodies by their domain
	std::vector<body_t> My_bodies(0); // Bodies that will be assigned to this processor
	
	double t_0 = (double)clock()/(double)CLOCKS_PER_SEC;
	ORB(nproc,bodies,bodies_in_orb_domains);
	CommunicateDomainsBodies(bodies,bodies_in_orb_domains,My_bodies);
	double t_1 = (double)clock()/(double)CLOCKS_PER_SEC;
	t_BodiesGeneration=t_BodiesGeneration+t_1-t_0;	//Add partitioning time to the time of bodies generation
	
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Build Tree (based on the global extent domain)
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−	
	domain_t G_domain;
	double G_dimmin,G_dimmax;
	// Get the extent over the bodies on all processors (global extent domain)
	getPointExtent(My_bodies,G_domain,G_dimmin,G_dimmax,true);
	// Build the tree by adding the local bodies
	TreeNode* root=Build_Tree(My_bodies,G_dimmin,G_dimmax);
	// Complete the local tree by adding the needed COMs (or bodies) on other processors
	CompleteLocalTrees(root,My_bodies,theta);
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Compute the forces on the bodies assigned to this processor 
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	double t_2 = (double)clock()/(double)CLOCKS_PER_SEC;
	std::vector<double3_t> force_BH(My_bodies.size());
	std::vector<double> force_Mag(My_bodies.size(),0);
	
	for ( int i = 0; i < My_bodies.size(); i++)
	{
		root->computeForceOnBody(My_bodies[i],theta,force_BH[i]);
		force_Mag[i] = sqrt((force_BH[i].r[0]*force_BH[i].r[0])+ 
							(force_BH[i].r[1]*force_BH[i].r[1]) + 
							(force_BH[i].r[2]*force_BH[i].r[2]));
		
	}
	double t_3 = (double)clock()/(double)CLOCKS_PER_SEC;

	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Display results 
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	if (rank==0)
	{
		
		std::cout<<"---------------------------------------"<<std::endl;
		std::cout<<"RESULTS (PARALLEL BARNES-HUT METHOD)"<<std::endl;
	}
	DisplayForceStat(force_Mag,t_BodiesGeneration,t_2-t_1,t_3-t_2,true);
		
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Verify the results
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	if (Verify==1)
	{
		std::vector<body_t> AllBodies(N*nproc); // Contians all the bodies generated on all processors
		// Collect all the bodies on rank 0 processor
		MPI_Gather(&bodies[0],N*sizeof(bodies[0]),MPI_BYTE,&AllBodies[0],
					   N*sizeof(bodies[0]),MPI_BYTE,0,MPI_COMM_WORLD);
		// Processor 0 calculate the forces using the serial implementation
		if (rank==0)
		{
			Serial_BH(AllBodies,theta,Verify,(t_1-t_0));
		}
	}	
	
	//−−−−−−−−−−−−−−−−−
	// Clean memory
	//−−−−−−−−−−−−−−−−−
	root->prune();
	root->~TreeNode();
	
}

//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Parse command line arguments
// N: Number of bodies of bodies per processor,
// Theta: Control parameter for Barnes-Hut algorithm
// Verify: Verify results using serial implementation of Barnes-Hut algorithim 
//			and using the direct method for force calculation
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
int interpret_arg(int argc, char** argv,int& N,double& theta,int& Verify) 
{  
	//−−−−−−−−−−−−−
    // MPI setup 
    //−−−−−−−−−−−−−
	int rank;	// Rank of the running processor
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);	//Identify the ID of the processor
  
	if (argc < 4)
	{
		if(rank==0)
		{	
			std::cerr << "You must pass :"<<std::endl
				<<"1.The number of bodies per processor"<<std::endl
				<<"2. Theta (Barnes-Hut algorithm control parmeter)"<<std::endl
				<<"3. Verify results (pass '1' to do verification)"
				<< std::endl;
		}
		
		MPI_Finalize();
		exit (1);
	}
	N = atoi(argv[1]) ;
	theta=strtod(argv[2],NULL);	
	Verify =atoi(argv[3]) ;	
}
	
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
// Generate 'N' random bodies and calculate the forces between them using Barnes-Hut algorithm
//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
void BarnesHutForceCalculation(int N, double theta, int Verify)
{
	//−−−−−−−−−−−−−
    // MPI setup 
    //−−−−−−−−−−−−−
	int rank, nproc;	// Rank of the running processor, Number of running processors
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);	//Get the number of processors available
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);	//Identify the ID of the processor
	
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Generate N random bodies on each processor
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	double t_0 = (double)clock()/(double)CLOCKS_PER_SEC;	// Start timing
	std::vector<body_t> bodies(N);
	Generate_Bodies(bodies);
	double t_1 = (double)clock()/(double)CLOCKS_PER_SEC;	// Record time for points generation
	
	if (nproc==1)
		Serial_BH(bodies,theta,Verify,(t_1-t_0));
	else
		Parallel_BH(bodies,N,theta,Verify,(t_1-t_0));

}

int main(int argc, char** argv)
{

    int rank, nproc;	// Rank of the running processor, Number of running processors
    int N,Verify;		// Number of bodies of bodies per processor, Verify results(1)
	double theta;		// Control parameter for Barnes-Hut algorithm

	//−−−−−−−−−−−−−
    // MPI setup 
    //−−−−−−−−−−−−−
	MPI_Init(&argc, &argv);					//Initialize MPI Environment
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);	//Get the number of processors available
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);	//Identify the ID of the processor
    
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	// Parse command line arguments
	//−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
	interpret_arg(argc,argv,N,theta,Verify);
	
	BarnesHutForceCalculation(N,theta,Verify);

	
	MPI_Finalize();
    return 0;
}

