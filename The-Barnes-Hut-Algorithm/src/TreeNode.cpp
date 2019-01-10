//-------------------------------------------------
// Parallel Implementation of the Barnes-Hut algorithm using OpenMPI Library
// Course: 'ECE4530-Parallel Processing' at the University of Manitoba.
// Implemented by: Micheal Sobhy
// (The initial code structure and some functions were provided by the course instructor)
// Email: sobhymich@gmail.com
//----------------------------------------------------

#include <iostream>
#include "TreeNode.h"
#include <math.h>
#include <cassert>

//----------------------------------------------------
//Construct a Tree Node, the domain is defined by the min and max point 
//----------------------------------------------------
TreeNode::TreeNode(double x1, double x2, double y1, double y2, double z1, double z2)
{
    xmin_ = x1;
    ymin_ = y1;
    zmin_ = z1;
    xmax_ = x2;
    ymax_ = y2;
    zmax_ = z2;
	
    //Initialize the tree Node (NULL Children & 0 Bodies)
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

//----------------------------------------------------
// Divide a Tree Node domain into 8 sub-domains by creating 8 children nodes
//----------------------------------------------------
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

//----------------------------------------------------
//Check if the body belong to this tree node (in the domain of the node)
//----------------------------------------------------
bool TreeNode::containsBody(const body_t& body) const
{
    if (body.r[0] >= xmin_ && body.r[0] <= xmax_ && body.r[1] >= ymin_ && body.r[1] <= ymax_ && body.r[2] >= zmin_ && body.r[2] <= zmax_) return true;
    return false;
}
//----------------------------------------------------
// Add a body to the tree such that every tree leaf contains only one body
// and the bodies are place in the nodes corresponding to thier domain
//----------------------------------------------------
void TreeNode::addBody(const body_t& body)
{
	//The node is a leaf & EMPTY (push the body)
    if (number_of_bodies_ == 0)
    {
        body_ = body;
        number_of_bodies_ = 1;
    }
	//The node is a leaf & CONTAIN a body (create children nodes and push the original and pushed body to the new leafs
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


void TreeNode::getCoM(body_t& com) const
{
    com = com_;
}
//----------------------------------------------------
//Calculate the Center of Mass (COM) for this tree node 
//----------------------------------------------------
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


//----------------------------------------------------
// Calculate the force applied on 'body' due to the bodies in this tree node
//----------------------------------------------------
void TreeNode::computeForceOnBody(const body_t& body, double theta, double3_t& F) const
{
	//----------------------------------------------------
	// Tree Node is not a leaf.
	// Check if we can use the COM to calculate the force on the body, 
	// otherwise calculate the forces using the children tree nodes
	//----------------------------------------------------
    if (number_of_bodies_ > 1)
    {
		//Calculate D: This tree node dimension (size of its domain)
		double D = sqrt((xmax_- xmin_)*(xmax_-xmin_)+(ymax_- ymin_)*(ymax_-ymin_)+(zmax_- zmin_)*(zmax_-zmin_));	
		//Calculate R: Distance between 'body' and the COM of this tree node
		double R = sqrt((body.r[0] - com_.r[0])*(body.r[0] - com_.r[0]) + (body.r[1] - com_.r[1])*(body.r[1] - com_.r[1]) + (body.r[2] - com_.r[2])*(body.r[2] - com_.r[2]));
		
		// Body is far enough from the domain of this tree node
		// (use th COM of that node to calculate the force)
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
		// Body is not far enough from the domain of this tree node
		// Calculate the forces using the children tree nodes
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
	
	//----------------------------------------------------
	// Tree node is a leaf and contains a body.
	// Calculate the force due to the body in the leaf node.
	//----------------------------------------------------
	else if (number_of_bodies_ == 1)
    {
        // note com is equal to body in this case by our convention.
        double R = sqrt((body.r[0] - com_.r[0])*(body.r[0] - com_.r[0]) + (body.r[1] - com_.r[1])*(body.r[1] - com_.r[1]) + (body.r[2] - com_.r[2])*(body.r[2] - com_.r[2]));
        
        if (R > 1e-14)
        {
            double tmp = G*body.m*com_.m/(R*R*R);
            F.r[0] -= tmp*(body.r[0] - com_.r[0]);
            F.r[1] -= tmp*(body.r[1] - com_.r[1]);
            F.r[2] -= tmp*(body.r[2] - com_.r[2]);
        }
    }
	//----------------------------------------------------
	// Tree node is a leaf and empty.
	// No force to be applied on the body.
	//----------------------------------------------------
    else return;
}


//----------------------------------------------------
// Each processor needs to have a complete tree, i.e., contains information about bodies 
// in all domains, to be able to calculate the force due to all the bodies.
// This function finds the locally essential COMs (or bodies) in my tree required
// to complete the tree on a processor computing the forces for the bodies in a given domain 
//----------------------------------------------------
void TreeNode::LETBodies(const domain_t& domain, double theta, std::vector<body_t>& bodies)
{
	int rank, nproc;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
	//----------------------------------------------------
	// Tree node is an empty leaf
	// No bodies
	//----------------------------------------------------
    if (number_of_bodies_ == 0) return;
	
	//----------------------------------------------------
    // Tree node is a leaf and contains a body.
	// Push the body in the list to be sent to other processor
	//----------------------------------------------------
	else if (number_of_bodies_ == 1)
    {
        bodies.push_back(body_);
        return;
    }
    
	//----------------------------------------------------
	// Tree node is not a leaf 
	// If the domain is far enough push the COM information,
	// otherwise, check the children tree nodes
	//----------------------------------------------------
    
    

    // Calculate the minimum distance between the current node and the given domain.
    double R2 = 0;
        
    if (domain.max[0] < xmin_) R2 += (domain.max[0] - xmin_)*(domain.max[0] - xmin_);
    else if (domain.min[0] > xmax_) R2 += (domain.min[0] - xmax_)*(domain.min[0] - xmax_);
        
    if (domain.max[1] < ymin_) R2 += (domain.max[1] - ymin_)*(domain.max[1] - ymin_);
    else if (domain.min[1] > ymax_) R2 += (domain.min[1] - ymax_)*(domain.min[1] - ymax_);
        
    if (domain.max[2] < zmin_) R2 += (domain.max[2] - zmin_)*(domain.max[2] - zmin_);
    else if (domain.min[2] > zmax_) R2 += (domain.min[2] - zmax_)*(domain.min[2] - zmax_);
        
    double R = sqrt(R2);

	//Calculate D: This tree node dimension (size of its domain)
	double D = sqrt((xmax_- xmin_)*(xmax_-xmin_)+(ymax_- ymin_)*(ymax_-ymin_)+(zmax_- zmin_)*(zmax_-zmin_));	
	
	if (R!=0)
	{	
		// Current tree node is far enough, push the COM in the list to be sent to other processor 
		if (D/R<=theta)
		{
			bodies.push_back(com_);
		}
		// Current tree node is not far enough, check the children tree nodes
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
	// Current tree node is so close, add the bodies in the children tree nodes
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



