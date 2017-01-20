// smallworld.cpp : Defines the entry point for the console application.
//

/*
This file uses part of ``kdtree'', a library for working with kd-trees:

Copyright (C) 2007-2009 John Tsiombikas <nuclear@siggraph.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/

#include <string.h>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>


using std::cerr;
using std::cout;
using std::string;
using std::vector;
using std::ifstream;
using std::exception;
using std::swap;
using std::lower_bound;
using std::nth_element;


/////////////////////////////////////////////////////////////////////////////////

enum { DIM = 2 };

inline double SQ(double x) { return x * x; }

struct kdnode
{
    double pos[DIM];
    int data;
	kdnode *left, *right;
};

typedef kdnode FriendInfo;


struct kdhyperrect
{
	/* minimum/maximum coords */
	double min[DIM], max[DIM];

	kdhyperrect(const double *min_, const double *max_)
	{
		memcpy(min, min_, sizeof(min));
		memcpy(max, max_, sizeof(max));
	}
};

/////////////////////////////////////////////////////////////////////////////////

inline bool hyperrect_dist_sq(kdhyperrect *rect, const double *pos, double limit)
{
	for (int i=0; i < DIM; i++) {
		if (pos[i] < rect->min[i]) {
			limit -= SQ(rect->min[i] - pos[i]);
		} else if (pos[i] > rect->max[i]) {
			limit -= SQ(rect->max[i] - pos[i]);
		}
        else
            continue;

        if (limit <= 0)
            return false;
	}

	return true;
}


void hyperrect_extend(kdhyperrect *rect, const double *pos)
{
	for (int i=0; i < DIM; i++) {
		if (pos[i] < rect->min[i]) {
			rect->min[i] = pos[i];
		}
		if (pos[i] > rect->max[i]) {
			rect->max[i] = pos[i];
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////


typedef vector<FriendInfo> FriendInfos;

template<int m_idx>
struct IsInfoLess
{
    bool operator() (const FriendInfo& x, const FriendInfo& y)
    {
        return x.pos[m_idx] < y.pos[m_idx];
    }
    bool operator() (const FriendInfo* x, const FriendInfo* y)
    {
        return x->pos[m_idx] < y->pos[m_idx];
    }
};

template <int dir>
kdnode* insert(vector<FriendInfo*>::iterator begin, vector<FriendInfo*>::iterator end)
{
    if (begin == end)
        return 0;

    int diff = int(end - begin);
    if (1 == diff)
    {
		kdnode* node = *begin;

        node->left = 0;
        node->right = 0;

        return node;
    }

    int halfSize = diff / 2;

    vector<FriendInfo*>::iterator middle = begin + halfSize;

	nth_element(begin, middle, end, IsInfoLess<dir>());

	kdnode* node = *middle;

	enum { new_dir = (dir + 1) % DIM };
    node->left = insert<new_dir>(begin, middle);
    node->right = insert<new_dir>(++middle, end);

    return node;
}


struct SearchResult
{
    double dist_sq;
    kdnode* node;
    SearchResult(double dist_sq_, kdnode* node_)
        : dist_sq(dist_sq_), node(node_) {}
};

bool operator < (const SearchResult& left, const SearchResult& right)
{
    return left.dist_sq < right.dist_sq;
}

bool operator < (const SearchResult& left, double right)
{
    return left.dist_sq < right;
}

bool operator < (double left, const SearchResult& right)
{
    return left < right.dist_sq;
}

template<int dir>
void kd_nearest_i(kdnode *node, const double *pos,
				  vector<SearchResult>& result, kdhyperrect* rect)
{
	kdnode *nearer_subtree, *farther_subtree;
	double *nearer_hyperrect_coord, *farther_hyperrect_coord;

	/* Decide whether to go left or right in the tree */
	double dummy = pos[dir] - node->pos[dir];
	if (dummy <= 0) {
		nearer_subtree = node->left;
		farther_subtree = node->right;
		nearer_hyperrect_coord = rect->max + dir;
		farther_hyperrect_coord = rect->min + dir;
	} else {
		nearer_subtree = node->right;
		farther_subtree = node->left;
		nearer_hyperrect_coord = rect->min + dir;
		farther_hyperrect_coord = rect->max + dir;
	}

	enum { new_dir = (dir + 1) % DIM };

	if (nearer_subtree) {
		/* Slice the hyperrect to get the hyperrect of the nearer subtree */
		dummy = *nearer_hyperrect_coord;
		*nearer_hyperrect_coord = node->pos[dir];
		/* Recurse down into nearer subtree */
		kd_nearest_i<new_dir>(nearer_subtree, pos, result, rect);
		/* Undo the slice */
		*nearer_hyperrect_coord = dummy;
	}

    bool do_farther_subtree = false;
    if (farther_subtree != 0)
    {
		/* Get the hyperrect of the farther subtree */
		dummy = *farther_hyperrect_coord;
		*farther_hyperrect_coord = node->pos[dir];

		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our minimum distances. */
        do_farther_subtree = (result.size() < 3 
            || hyperrect_dist_sq(rect, pos, result.rbegin()->dist_sq));
    }

    if (node->pos != pos && (0 == farther_subtree || do_farther_subtree))
    {
	    /* Check the distance of the point at the current node, compare it with our bests so far */
	    double dist_sq = SQ(node->pos[0] - pos[0]);
	    for(int i = 1; i < DIM; i++) {
		    dist_sq += SQ(node->pos[i] - pos[i]);
	    }

        do
        {
            if (result.size() >= 3)
            {
                if (result.rbegin()->dist_sq < dist_sq)
                    break;
                result.pop_back();
            }
            vector<SearchResult>::iterator it = lower_bound(result.begin(), result.end(), dist_sq);
            result.insert(it, SearchResult(dist_sq, node));
        }
        while (false);
    }

	if (do_farther_subtree)
    {
		/* Recurse down into farther subtree */
		kd_nearest_i<new_dir>(farther_subtree, pos, result, rect);
	}

	if (farther_subtree != 0)
    {
		/* Undo the slice on the hyperrect */
		*farther_hyperrect_coord = dummy;
    }
}


template<int dir>
void kd_nearest_i_nearer_subtree(kdnode *node, const double *pos,
				  vector<SearchResult>& result, kdhyperrect* rect, int* flags)
{
	kdnode *nearer_subtree, *farther_subtree;
	double *nearer_hyperrect_coord, *farther_hyperrect_coord;
    int flag;

	/* Decide whether to go left or right in the tree */
	double dist = pos[dir] - node->pos[dir];
	if (dist <= 0) {
		nearer_subtree = node->left;
		farther_subtree = node->right;
		nearer_hyperrect_coord = rect->max + dir;
		farther_hyperrect_coord = rect->min + dir;
        flag = 1 << (dir * 2);
	} else {
		nearer_subtree = node->right;
		farther_subtree = node->left;
		nearer_hyperrect_coord = rect->min + dir;
		farther_hyperrect_coord = rect->max + dir;
        flag = 1 << (dir * 2 + 1);
	}

	enum { new_dir = (dir + 1) % DIM };

	if (nearer_subtree) {
		/* Slice the hyperrect to get the hyperrect of the nearer subtree */
		double dummy = *nearer_hyperrect_coord;
		*nearer_hyperrect_coord = node->pos[dir];
		/* Recurse down into nearer subtree */
		kd_nearest_i_nearer_subtree<new_dir>(nearer_subtree, pos, result, rect, flags);
		/* Undo the slice */
		*nearer_hyperrect_coord = dummy;
	}

    if (*flags & flag)
        return;

    bool do_farther_subtree = false;
    double dummy = 0;
    if (farther_subtree != 0)
    {
		/* Get the hyperrect of the farther subtree */
		dummy = *farther_hyperrect_coord;
		*farther_hyperrect_coord = node->pos[dir];

		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our minimum distances. */
        do_farther_subtree = (result.size() < 3 
            || SQ(dist) < result.rbegin()->dist_sq);

        if (!do_farther_subtree)
            *flags |= flag;
    }

    if (node->pos != pos && (0 == farther_subtree || do_farther_subtree))
    {
	    /* Check the distance of the point at the current node, compare it with our bests so far */
	    double dist_sq = SQ(node->pos[0] - pos[0]);
	    for(int i = 1; i < DIM; i++) {
		    dist_sq += SQ(node->pos[i] - pos[i]);
	    }

        do
        {
            if (result.size() >= 3)
            {
                if (result.rbegin()->dist_sq < dist_sq)
                    break;
                result.pop_back();
            }
            vector<SearchResult>::iterator it = lower_bound(result.begin(), result.end(), dist_sq);
            result.insert(it, SearchResult(dist_sq, node));
        }
        while (false);
    }

	if (do_farther_subtree)
    {
		/* Recurse down into farther subtree */
		kd_nearest_i<new_dir>(farther_subtree, pos, result, rect);
	}

	if (farther_subtree != 0)
    {
		/* Undo the slice on the hyperrect */
		*farther_hyperrect_coord = dummy;
    }
}

/////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		cerr << "Usage: smallworld input_file\n";
		return 1;
	}

	try
	{
	    //*
		string path(argv[1]);
		if (string::npos == path.find_first_of("\\/"))
		{
			string dir(argv[0]);
			string::size_type pos = dir.find_last_of("\\/");
			if (string::npos != pos)
			{
				path = dir.substr(0, pos + 1) + path;
			}
		}

		ifstream inFile(path.c_str());
		//*/
		//ifstream inFile(argv[1]);

        FriendInfos infos;

		while (!inFile.eof())
        {
            FriendInfo info;
            inFile >> info.data >> info.pos[0] >> info.pos[1];

			if (inFile.fail())
				break;

            infos.push_back(info);
        }

        if (infos.size() < 2)
            return 0;

        vector<FriendInfo*> infoPtrs;

        for (FriendInfos::iterator it = infos.begin(); it != infos.end(); ++it)
        {
            infoPtrs.push_back(&*it);
        }

        kdnode* root = insert<0>(infoPtrs.begin(), infoPtrs.end());

		kdhyperrect rect(infos.begin()->pos, infos.begin()->pos);

        for (FriendInfos::iterator it = infos.begin(); ++it != infos.end(); )
        {
            hyperrect_extend(&rect, it->pos);
        }

        vector<SearchResult> result;
		kdhyperrect tempRect = rect;
        // search
        for (FriendInfos::iterator it = infos.begin(); it != infos.end(); ++it)
        {
            result.clear();

			//kdhyperrect tempRect = rect;

            int flags = 0;

			/* Search for the nearest neighbours recursively */
			kd_nearest_i_nearer_subtree<0>(root, it->pos, result, &tempRect, &flags);

            /*
            vector<SearchResult>::iterator itRes = result.begin();
			cout << it->data << ' ' << itRes->node->data;
            while (++itRes != result.end())
            {
                cout << ',' << itRes->node->data;
            }
            cout << '\n';
            */
            printf("%d %d,%d,%d\n", it->data, result[0].node->data, result[1].node->data, result[2].node->data);
        }
	}
	catch (exception& e)
	{
		cerr << e.what();
		return 1;
	}

	return 0;
}

