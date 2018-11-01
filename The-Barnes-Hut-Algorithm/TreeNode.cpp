#include <iostream>
#include "TreeNode.h"
#include <math.h>
#include <cassert>

//Construct a Tree Node 
//defined by point min and point max 
//Initialize the tree Node (NULL Children & 0 Bodies)
TreeNode::TreeNode(double x1, double x2, double y1, double y2, double z1, double z2)
{
    xmin_ = x1;
    ymin_ = y1;
    zmin_ = z1;
    xmax_ = x2;
    ymax_ = y2;
    zmax_ = z2;
    
    number_of_bodies_ = 0;
    for (unsigned int ichild = 0; ichild < 8; ichild++) children_[ichild] = NULL;
}

TreeNode::~TreeNode()
{
    for (unsigned int ichild = 0; ichild < 8; ichild++)
    {
        if (children_[ichild] != NULL) delete children_[ichild];
        children_[ichild] = NULL;
    }
}
//Divde a box to 8 boxes 
//Set the Children nodes for a given tree node
void TreeNode::spawnChildren()
{
    double dx = 0.5*(xmax_ - xmin_);
    double dy = 0.5*(ymax_ - ymin_);
    double dz = 0.5*(zmax_ - zmin_);
    
    for (unsigned int ichild = 0; ichild < 8; ichild++)
    {
        assert(children_[ichild] == NULL);
    }
    
    children_[0] = new TreeNode(xmin_, xmin_ + dx, ymin_, ymin_ + dy, zmin_, zmin_ + dz);
    children_[1] = new TreeNode(xmin_ + dx, xmax_, ymin_, ymin_ + dy, zmin_, zmin_ + dz);
    children_[2] = new TreeNode(xmin_, xmin_ + dx, ymin_ + dy, ymax_, zmin_, zmin_ + dz);
    children_[3] = new TreeNode(xmin_ + dx, xmax_, ymin_ + dy, ymax_, zmin_, zmin_ + dz);
    children_[4] = new TreeNode(xmin_, xmin_ + dx, ymin_, ymin_ + dy, zmin_ + dz, zmax_);
    children_[5] = new TreeNode(xmin_ + dx, xmax_, ymin_, ymin_ + dy, zmin_ + dz, zmax_);
    children_[6] = new TreeNode(xmin_, xmin_ + dx, ymin_ + dy, ymax_, zmin_ + dz, zmax_);
    children_[7] = new TreeNode(xmin_ + dx, xmax_, ymin_ + dy, ymax_, zmin_ + dz, zmax_);
}
//Check if the body belong to this tree node (in the box of the node)
bool TreeNode::containsBody(const body_t& body) const
{
    if (body.r[0] >= xmin_ && body.r[0] <= xmax_ && body.r[1] >= ymin_ && body.r[1] <= ymax_ && body.r[2] >= zmin_ && body.r[2] <= zmax_) return true;
    return false;
}

//Function to Add a body to the tree where every tree node contains only one body 
//Reorder the tree where the bodies are on the leafs and placed in the correct node (box)
void TreeNode::addBody(const body_t& body)
{
	//The node is a leaf & EMPTY (push the body)
    if (number_of_bodies_ == 0)
    {
        body_ = body;
        number_of_bodies_ = 1;
    }
	//The node is a leaf & CONTAIN a body (create children nodes and push the original and pused body to the new leafs
    else if (number_of_bodies_ == 1)
    {
        spawnChildren();

        bool found_child = false;
        for (unsigned int ichild = 0; ichild < 8; ichild++)
        {
            if (children_[ichild]->containsBody(body_))
            {
                children_[ichild]->addBody(body_);
                ichild = 8; //break loop
                found_child = true;
            }
        }
        if (!found_child)
        {
            std::cerr << "Could not find a child for body ( " << body_.r[0] << ", " << body_.r[1] << ", " << body_.r[2] << " ) " << std::endl;
            std::cerr << "In box ( " << xmin_ << ", " << ymin_ << ", " << zmin_ << " ) to ( " << xmax_ << ", " << ymax_ << ", " << zmax_ << " ) " << std::endl;
        }

        found_child = false;
        for (unsigned int ichild = 0; ichild < 8; ichild++)
        {
            if (children_[ichild]->containsBody(body))
            {
                children_[ichild]->addBody(body);
                ichild = 8; //break loop
                found_child = true;
            }
        }
        if (!found_child)
        {
            std::cerr << "Could not find a child for body ( " << body_.r[0] << ", " << body_.r[1] << ", " << body_.r[2] << " ) " << std::endl;
            std::cerr << "In box ( " << xmin_ << ", " << ymin_ << ", " << zmin_ << " ) to ( " << xmax_ << ", " << ymax_ << ", " << zmax_ << " ) " << std::endl;
        }
        number_of_bodies_++;
    }
	//Node is NOT A LEAF (push the body to the correct node)
    else
    {
        bool found_child = false;
        for (unsigned int ichild = 0; ichild < 8; ichild++)
        {
            if (children_[ichild]->containsBody(body))
            {
                children_[ichild]->addBody(body);
                ichild = 8;
                found_child = true;
            }
        }
        number_of_bodies_++;
        
        if (!found_child)
        {
            std::cerr << "Could not find a child for body ( " << body_.r[0] << ", " << body_.r[1] << ", " << body_.r[2] << " ) " << std::endl;
            std::cerr << "In box ( " << xmin_ << ", " << ymin_ << ", " << zmin_ << " ) to ( " << xmax_ << ", " << ymax_ << ", " << zmax_ << " ) " << std::endl;
        }

    }
}

void TreeNode::prune()
{
    for (unsigned int ichild = 0; ichild < 8; ichild++)
    {
        if (children_[ichild] != NULL)
        {
            if (children_[ichild] -> getBodyCount() < 1)
            {
                delete children_[ichild];
                children_[ichild] = NULL;
            }
            else
            {
                children_[ichild]->prune();
            }
        }
    }
}

//starting from TREE LEVEL='level' GET maxlevel, nbodies, nnodes (below LEVEL) 
void TreeNode::diagnostics(int level, int& maxlevel, int& nbodies, int& nnodes) const
{
    if (level > maxlevel) maxlevel=level;
    nnodes++;
	
    bool isleaf = true;
    for (unsigned int ichild = 0; ichild < 8; ichild++)
    {
        if (children_[ichild] != NULL)
        {
            children_[ichild]->diagnostics(level+1, maxlevel, nbodies, nnodes);
            isleaf = false;
        }
    }
    if (isleaf)
    {
        nbodies+= number_of_bodies_;
    }
}

void TreeNode::getCoM(body_t& com) const
{
    com = com_;
}

//Calculate the COM for the tree node 
void TreeNode::computeCoM()
{
	//LEAF AND EMPTY (COM.m=0 COM.r = center of the square)
    if (number_of_bodies_ == 0)
    {
        com_.r[0] = 0.5*(xmax_ + xmin_);
        com_.r[1] = 0.5*(ymax_ + ymin_);
        com_.r[2] = 0.5*(zmax_ + zmin_);
        com_.m = 0.0;
    }
	//LEAF AND CONTAINS A BODY (COM = BODY)
    else if (number_of_bodies_ == 1)
    {
        com_ = body_;
    }
	//NOT A LEAF
    else
    {
		com_.r[0] = 0.0;
		com_.r[1] = 0.0;
		com_.r[2] = 0.0;
		com_.m = 0.0;
		body_t temp;
		bool all_null = true;
		for ( int ichild = 0; ichild < 8; ichild ++)
		{
			if (children_[ichild] !=NULL)
			{
				all_null=false;
				children_[ichild]->computeCoM();
				children_[ichild]->getCoM(temp);
				
				com_.m +=temp.m;
				com_.r[0] +=temp.r[0]*temp.m;
				com_.r[1] +=temp.r[1]*temp.m;
				com_.r[2] +=temp.r[2]*temp.m;
			}
		}
		
		assert(all_null ==false);
		com_.r[0] /=com_.m;
		com_.r[1] /=com_.m;
		com_.r[2] /=com_.m;
    }
}


void TreeNode::computeForceOnBody(const body_t& body, double theta, double3_t& F) const
{
	//NOT A LEAF
    if (number_of_bodies_ > 1)
    {
		//Calculate D
		double D = sqrt((xmax_- xmin_)*(xmax_-xmin_)+(ymax_- ymin_)*(ymax_-ymin_)+(zmax_- zmin_)*(zmax_-zmin_));	
		//Calculate R
		double R = sqrt((body.r[0] - com_.r[0])*(body.r[0] - com_.r[0]) + (body.r[1] - com_.r[1])*(body.r[1] - com_.r[1]) + (body.r[2] - com_.r[2])*(body.r[2] - com_.r[2]));
		
		//Point is Far enough (use th COM of that node to calculate the force)
		if (D/R<=theta)
		{
			if (R > 1e-14)
        	{
            	double tmp = G*body.m*com_.m/(R*R*R);
            	F.r[0] -= tmp*(body.r[0] - com_.r[0]);
            	F.r[1] -= tmp*(body.r[1] - com_.r[1]);
            	F.r[2] -= tmp*(body.r[2] - com_.r[2]);
        	}
		}
		//Calculate Forces using the children
		else
		{
			bool all_null = true;
			for ( int ichild = 0; ichild < 8; ichild ++)
			{
				if (children_[ichild] !=NULL)
				{
					all_null=false;
					children_[ichild]->computeForceOnBody(body,theta,F);
				}
			}
			assert(all_null ==false);
		}
 
    }
	//LEAF AND CONTAINS A BODY
    else if (number_of_bodies_ == 1)
    {
        //note com is equal to body in this case by our convention.
        double R = sqrt((body.r[0] - com_.r[0])*(body.r[0] - com_.r[0]) + (body.r[1] - com_.r[1])*(body.r[1] - com_.r[1]) + (body.r[2] - com_.r[2])*(body.r[2] - com_.r[2]));
        
        if (R > 1e-14)
        {
            double tmp = G*body.m*com_.m/(R*R*R);
            F.r[0] -= tmp*(body.r[0] - com_.r[0]);
            F.r[1] -= tmp*(body.r[1] - com_.r[1]);
            F.r[2] -= tmp*(body.r[2] - com_.r[2]);
        }
    }
	//LEAF AND EMPTY
    else return;
}

void TreeNode::LETBodies(const domain_t& domain, double theta, std::vector<body_t>& bodies)
{
	int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
	//LEAF AND EMPTY
    if (number_of_bodies_ == 0) return;
    //LEAF AND CONTAINS A BODY (PUSH IT)
	else if (number_of_bodies_ == 1)
    {
        bodies.push_back(body_);
        return;
    }
    
	//NOT A LEAF 
    //no need for an else here as the above conditions return.
    
    //------------------------------------------------------------------------
    //The following code will compute the minimum distance between the current
    //node and the given domain.
    //------------------------------------------------------------------------
    
    double R2 = 0;
        
    if (domain.max[0] < xmin_) R2 += (domain.max[0] - xmin_)*(domain.max[0] - xmin_);
    else if (domain.min[0] > xmax_) R2 += (domain.min[0] - xmax_)*(domain.min[0] - xmax_);
        
    if (domain.max[1] < ymin_) R2 += (domain.max[1] - ymin_)*(domain.max[1] - ymin_);
    else if (domain.min[1] > ymax_) R2 += (domain.min[1] - ymax_)*(domain.min[1] - ymax_);
        
    if (domain.max[2] < zmin_) R2 += (domain.max[2] - zmin_)*(domain.max[2] - zmin_);
    else if (domain.min[2] > zmax_) R2 += (domain.min[2] - zmax_)*(domain.min[2] - zmax_);
        
    double R = sqrt(R2);
/*	
	std::cout<<"I'm rank "<<rank<<" R = "<<R<<std::endl;
	std::cout<<" domain.max[0] "<<domain.max[0]<<" domain.max[1] "<<domain.max[1]<<" domain.max[2] "<<domain.max[2]<<std::endl;
	std::cout<<" xmin_ "<<xmin_<<" ymin_ "<<ymin_<<" zmin_ "<<zmin_<<std::endl;
	std::cout<<" xmax_ "<<xmax_<<" ymax_ "<<ymax_<<" zmax_ "<<zmax_<<std::endl;
*/
	//Calculate D
	double D = sqrt((xmax_- xmin_)*(xmax_-xmin_)+(ymax_- ymin_)*(ymax_-ymin_)+(zmax_- zmin_)*(zmax_-zmin_));	
	
	if (R!=0)
	{	
		//Current TREE NODE is Far enough (use it's COM)
		if (D/R<=theta)
		{
			bodies.push_back(com_);
		}
		//Not Far enough check the children
		else
		{
			bool all_null = true;
				for ( int ichild = 0; ichild < 8; ichild ++)
				{
					if (children_[ichild] !=NULL)
					{
						all_null=false;
						children_[ichild]->LETBodies(domain, theta, bodies);
					}
				}
				assert(all_null ==false);
		}	
	}
	else
	{
		bool all_null = true;
		for ( int ichild = 0; ichild < 8; ichild ++)
		{
			if (children_[ichild] !=NULL)
			{
				all_null=false;
				children_[ichild]->LETBodies(domain, theta, bodies);
			}
		}
		assert(all_null ==false);
		
	}
}



