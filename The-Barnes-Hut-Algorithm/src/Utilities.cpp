//-------------------------------------------------
// Parallel Implementation of the Barnes-Hut algorithm using OpenMPI Library
// Course: 'ECE4530-Parallel Processing' at the University of Manitoba.
// Implemented by: Micheal Sobhy
// (The initial code structure and some functions were provided by the course instructor)
// Email: sobhymich@gmail.com
//----------------------------------------------------

#include <vector>
#include <math.h>
#include <cassert>
#include <mpi.h>
#include <queue>
#include "./TreeNode.h"
#include <algorithm>

using namespace std;

//-----------------------------------------------------
// For a list of bodies, compute the domain that
// contains them (tight bound) in domain.
// Also return the minimum and maximum location over all
// three dimensions in dimmin and dimmax. The global flag
// can be used to make the domains local (false) or
// global (true).
//-----------------------------------------------------
void getPointExtent(const std::vector<body_t>& bodies, domain_t& domain, double& dimmin, double& dimmax, bool global)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    double local_xmin = HUGE_VAL;
    double local_xmax = -HUGE_VAL;
    double local_ymin = HUGE_VAL;
    double local_ymax = -HUGE_VAL;
    double local_zmin = HUGE_VAL;
    double local_zmax = -HUGE_VAL;
    
    for (unsigned int ibody = 0; ibody < bodies.size(); ibody++)
    {
        if (bodies[ibody].r[0] < local_xmin) local_xmin = bodies[ibody].r[0];
        if (bodies[ibody].r[0] > local_xmax) local_xmax = bodies[ibody].r[0];
        if (bodies[ibody].r[1] < local_ymin) local_ymin = bodies[ibody].r[1];
        if (bodies[ibody].r[1] > local_ymax) local_ymax = bodies[ibody].r[1];
        if (bodies[ibody].r[2] < local_zmin) local_zmin = bodies[ibody].r[2];
        if (bodies[ibody].r[2] > local_zmax) local_zmax = bodies[ibody].r[2];
    }
    
    double local_dimmin = min(min(local_xmin,local_ymin),local_zmin);
    double local_dimmax = max(max(local_xmax,local_ymax),local_zmax);
    
    double xmin, xmax, ymin, ymax, zmin, zmax;
    
    if (global == true)
    {
        MPI_Allreduce(&local_xmin, &xmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_ymin, &ymin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_zmin, &zmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        
        MPI_Allreduce(&local_xmax, &xmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&local_ymax, &ymax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(&local_zmax, &zmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        
        MPI_Allreduce(&local_dimmin, &dimmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        MPI_Allreduce(&local_dimmax, &dimmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    else
    {
        xmin = local_xmin;
        xmax = local_xmax;
        ymin = local_ymin;
        ymax = local_ymax;
        zmin = local_zmin;
        zmax = local_zmax;
        dimmin = local_dimmin;
        dimmax = local_dimmax;
    }
    
    domain.min[0] = xmin;
    domain.min[1] = ymin;
    domain.min[2] = zmin;
    
    domain.max[0] = xmax;
    domain.max[1] = ymax;
    domain.max[2] = zmax;
    
}


//-----------------------------------------------------
// Extract a single coordinate list from a body_t list
//-----------------------------------------------------
std::vector<double> getSingleCoordinateListFromBodies(const std::vector<body_t>& bodies, int dim)
{
    std::vector<double> coordinate_list(bodies.size());
    if (dim > 3 || dim < 0)
    {
        std::cerr << "Requested dimension " << dim << " is out of bounds at line " << __LINE__ << " of file " << __FILE__ << " for function " << __FUNCTION__ << std::endl;
        return coordinate_list;
    }
    
    for (unsigned int ipoint = 0; ipoint < bodies.size(); ipoint++)
    {
        coordinate_list[ipoint] = bodies[ipoint].r[dim];
    }
    
}


//-----------------------------------------------------
// For a list of bodies, partition them into 'P' orthogonal 
// domains using Orthogonal Recursive Bisection (ORB)
// partitioning method.
// Return the paritioned bodies in 'bodies_in_orb_domains' vector of vectors,
// where bodies_in_orb_domains[iproc] stores the local element
// indeces of all processors in subdomain iproc (i.e., destined for iproc).
//-----------------------------------------------------
void ORB(int P, const std::vector<body_t>& bodies, std::vector<std::vector<int> >& bodies_in_orb_domains)
{
	
    int rank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	//Queue of Domains to be partittioned 
    std::queue<Domain> Domains;
	//DI: initial domain (contain the bodies first read on the processor)
	//D0: Domain to be partitioned 
	//D_L: Domain of bodies before the weighted median
	//D_R: Domain of bodies after the weighted medain
	Domain DI,D0,D_L,D_R;
	//list and sorted list of a certain Coordinate for the bodies in the domain
	std::vector <double> list,sorted_list;
	
	//W_L,W_R new weights for the left and right domain
	//W_median: position of the wighted median with respect to the bodies on all domains
	//W_median_proc, W_median_proc_index: the processor where the W_median exist and it's index
	//dim: Coordinate to be used for division
	//SL_size: size of the sorted list on each processor
	//cutoff: index of the first element with respect to the global sorted list
	int W_L,W_R,W_median,W_median_proc,W_median_proc_index,dim,SL_size,cutoff;
	//for checking
	int ec1;
	//The margin to be used for division 
	double W_median_val;
	//All_SL_sizes: Array for the values of SL_size on all processor
	//cuttoffs: Array for the values of cutoff on all processors (last element in the array is the sum of number of elements on all processors)
	int All_SL_sizes[nproc],cutoffs[nproc+1];
	
	//Set ids for the bodies in the initail domain
	DI.pids.resize(bodies.size());
	for(int i=0;i<bodies.size();i++) {DI.pids[i]=i;}	
	
	//Set initial coordinate, weight, and bodies
	DI.Dim=0;	//Dims	0 1 2	
	DI.W=P;
	DI.bodies=bodies;
	//Push the initail domain to the queue
	Domains.push(DI);
	
	//if there is more domains in the queue to be partitioned
	while (!Domains.empty())
	{
		
		//Set my current domain from the queue
		D0=Domains.front();
		Domains.pop();
		dim=D0.Dim+1;
		//return to first dim (0) if we passed dim(2)
		if (dim>2) {dim=0;}
		
		//Get the list of bodies based on dim and sort them
		list.resize(0);
		sorted_list.resize(0);
		list=getSingleCoordinateListFromBodies(D0.bodies,dim);
		parallelBucketSort(list,sorted_list);
		
		//Calculate new weights (L,R)
		W_L=floor(D0.W/2);
		W_R=ceil(D0.W/2);

		//I will make processor 0 get the median element
		//Get the median elementlist
		//get the size of the sorted list on each processor
		SL_size=sorted_list.size();
		//Gather the SL_size from all processors to All_SL_sizes on rank 0
		MPI_Gather(&SL_size,1,MPI_INT,&All_SL_sizes[0],1,MPI_INT,0,MPI_COMM_WORLD);
	
		//---------------Rank 0------------------//
		if (rank==0)
		{
			
			W_median_proc=-1;
			ec1=0;
			cutoff=0;
			
			//Set cuttoffs array
			for(int i=0;i<nproc;i++)
			{
				cutoffs[i]=cutoff;
				cutoff+=All_SL_sizes[i];
			}
			cutoffs[nproc]=cutoff;
			
			//Calculate the weighted median 
			W_median=ceil((double)cutoffs[nproc]*(double)(W_L)/((double)W_L+(double)W_R));
			
			//Get the processor number where the W_median exist
			for(int i=0;i<nproc;i++)
			{
				if (W_median>=cutoffs[i] && W_median<cutoffs[i+1]) 
				{
					W_median_proc=i;
					ec1++;
				}
					
			}
			assert(ec1==1 && W_median_proc<nproc);
			//Get the index of the weighted median on W_median_proc
			W_median_proc_index=W_median-cutoffs[W_median_proc];
			assert(W_median_proc_index<=All_SL_sizes[W_median_proc]);
			
			//Communicate the value of the W_median_proc, and W_median_proc_index to all processors
			MPI_Bcast(&W_median_proc,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(&W_median_proc_index,1,MPI_INT,0,MPI_COMM_WORLD);
			
			//If rank 0 is W_median_proc get W_median_val
			if (W_median_proc==rank)
			{
				W_median_val=sorted_list[W_median_proc_index];
			}
			//Receive W_median_val from any processor 
			else
			{
				MPI_Recv(&W_median_val,1,MPI_DOUBLE,W_median_proc,100,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
			
			//Communicate W_median_val to all processors
			MPI_Bcast(&W_median_val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
		}
		//------------ALL PROCESSORS EXCEPT RANK 0----------//
		else
		{
			//Receive W_median_proc, and W_median_proc_index from rank 0
			MPI_Bcast(&W_median_proc,1,MPI_INT,0,MPI_COMM_WORLD);
			MPI_Bcast(&W_median_proc_index,1,MPI_INT,0,MPI_COMM_WORLD);
			
			//If W_median is on this rank get it and send it to rank 0
			if (W_median_proc==rank)
			{
				MPI_Send(&sorted_list[W_median_proc_index],1,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
			}
			
			//Receive W_median_val from rank 0
			MPI_Bcast(&W_median_val,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		}
	
	//---------FOR ALL PROCESSORS-----------//
	//Now all processors Hit this point with the value of the weighted median 
	//I need to divide the bodies to two domains then push the point to the queue
		
	//if R is true push the element on the W_median to the right (First W_median belong to the right domain)
	bool R=true;
		
		
	for (int i=0;i<D0.bodies.size();i++)
	{
		//Elements on the W_median push one to the left and one to the right 
		if (list[i]==W_median_val)
		{
			if (R)
			{
				D_R.bodies.push_back(D0.bodies[i]);
				D_R.pids.push_back(D0.pids[i]);
			}
			else
			{
				D_L.bodies.push_back(D0.bodies[i]);
				D_L.pids.push_back(D0.pids[i]);
			}
			R=!R;
		
		}
		
		//Put bodies before the W_median to the left domain
		else if (list[i]<W_median_val)
		{
			D_L.bodies.push_back(D0.bodies[i]);
			D_L.pids.push_back(D0.pids[i]);
		}
		//Put bodies after the W_median to the right domain
		else
		{
			D_R.bodies.push_back(D0.bodies[i]);
			D_R.pids.push_back(D0.pids[i]);
		}
	}
	
	//Set the parameters of the new domains
	D_L.W=W_L;
	D_L.Dim=dim;
	D_R.Dim=dim;
	D_R.W=W_R;
	
	//If the weight of the right domain is more than one push to queue in order to be repartitioned 
	if (W_L>1)
	{
		Domains.push(D_L);
	}
	// no more partitioning required store the ids for this domain in bodies_in_orb_domains
	else
	{
		bodies_in_orb_domains.push_back(D_L.pids);
	}
	//If the weight of the right domain is more than one push to queue in order to be repartitioned 
	if (W_R>1)
	{
		Domains.push(D_R);
	}
	// no more partitioning required store the ids for this domain in bodies_in_orb_domains
	else
	{
		bodies_in_orb_domains.push_back(D_R.pids);
	}
	
	//Resize to zero to be reused
	D_L.bodies.resize(0);
	D_L.pids.resize(0);
	D_R.bodies.resize(0);
	D_R.pids.resize(0);
		
	}
	
}


//-----------------------------------------------------
// Parallel approach for Bucket Sort.
//
// Inputs: The values to sort on each processor
//
// Outputs: A subset of the sorted values. 
//
// Note that the entries
// in sorted values are not the same as values to sort
// on any given processor. They are a subset of the
// total set of sorted values where rank 0 will contain
// the lowest sorted numbers.
//-----------------------------------------------------
void parallelBucketSort(const std::vector<double>& values_to_sort, std::vector<double>& sorted_values)
{
    
    int rank, nproc;
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
	
	double local_min,local_max,global_min,global_max;
	//P small buckets on each processor
	std::vector<std::vector<double> > S_Buckets(nproc);
	//Buckets from all processors to form the big  bucket
	std::vector<std::vector<double> > B_Buckets(nproc);
	//The Big sorted bucket on this processor
	std::vector<double> My_Bucket(0);
	//Value from 0 -> nproc (to decide the bucket number)
	double B_num;
	
	//Get local minimum and maximum
	local_min= *std::min_element(values_to_sort.begin(),values_to_sort.end());
	local_max= *std::max_element(values_to_sort.begin(),values_to_sort.end());
	//Get global minimum and maximum
	MPI_Allreduce(&local_min,&global_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);	
	MPI_Allreduce(&local_max,&global_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	
	
	for (int i=0;i<values_to_sort.size();i++)
	{
		//Calculate B_num
		B_num=(values_to_sort[i]-global_min)*nproc/(global_max-global_min);
		
		// Get the bucket number
		if (B_num<nproc)
		{
			B_num=floor(B_num);
		}
		//handle case when B_num== nproc
		else
		{
			B_num=nproc-1;
		}
		//Push the element to correct bucket
		S_Buckets[B_num].push_back(values_to_sort[i]);
	}
	
	//Communicate all buckets to get the buckets that belong to each processor
	MPI_Alltoall_vecvecT(S_Buckets,B_Buckets);

	//Put all the elements in 1 bucket 'My_Bucket'
	for(int k=0;k<nproc;k++)
	{
		for(int j=0;j<B_Buckets[k].size();j++)
		{
			My_Bucket.push_back(B_Buckets[k][j]);
		}
	}
	
	//Sort the elements in My_bucket
	sort(My_Bucket.begin(),My_Bucket.end());
	sorted_values=My_Bucket;
	
		
}


