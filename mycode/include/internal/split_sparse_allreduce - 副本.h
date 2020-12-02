#ifndef SPARSE_SPLIT_SCANNTER_ALLREDUCE_H
#define SPARSE_SPLIT_SCANNTER_ALLREDUCE_H

#include <map>
#include <vector>
#include <cstring>
#include<cmath>
#include<iostream>

#include "logging/simple_logging.h"
#include "internal/common.h"
#include "internal/reduce_operator.h"
#include "internal/simple_allreduce.h"
#include "internal/p2p_communication.h"


namespace spar {
//hierarchiar split-allgather的算法实现,allgather阶段先节点外通信再节点间通信
template<class O, class T>
void HierSplitAlgather(T *buffer, int count, int id, std::vector<int> &worker_list, MPI_Comm comm) {
int worker_number = worker_list.size();
    if (worker_number == 1) return;
    if (count < 512 * worker_number) {
        SimpleAllreduce<O>(buffer, count, id, worker_list, comm);
        return;
    }

    
	int block_size=count/worker_number;
	std::map<int,std::pair<int,int>> blocks; //记录每块的起始地址和块大小
	for(int i=0;i<worker_number;i++)
    {
		blocks[i].first=i*block_size;
	    blocks[i].second=block_size;
	}		
	blocks[worker_number-1].second = count - (worker_number - 1) * block_size;
	int right=-1;
	int left=-1;
	int my_index=-1;
	
	for(int i=0;i<worker_number;i++)
	{
		if(worker_list[i]==id)
		{
			my_index=i;
			break;
		}
	}
	 
	bool need_check = true;
	int isize = sizeof(int);
    int vsize = sizeof(T);
	T *send_buffer=new T[count];
	T *value_buffer = new T[count];
	int *index_buffer=new int[count];
	 
	T *recv_buffer=new T[count];
	//MPI_Request req[2];
	//MPI_Status status[2];
 
	//split 阶段
    for(int i=0;i<worker_number-1;++i)
	{
		left=worker_list[(my_index-i-1+worker_number)%worker_number];
	 	right=worker_list[(my_index+i+1)%worker_number];
		MPI_Status statuses[2];
        MPI_Request requests[2];
		T *base=buffer+blocks[right].first;
		block_size=blocks[right].second;
 
		if(need_check)//是否过滤标志
		{
			int nnz=0;
			for(int j=0;j<block_size;j++)
			{
				if(base[j]!=0)
				{
					++nnz;
				}
			}
			 
			if(nnz*(isize+vsize)<block_size * vsize)//满足稀疏传送条件
			{
				int k=0;
				for(int j=0;j<block_size;++j)
				{
					if(base[j]!=0)
					{
						value_buffer[k]=base[j];
						index_buffer[k]=j;
						++k;
					}					
				}
				memcpy(value_buffer+nnz,index_buffer,isize*nnz); //将对应索引值放在value值的后面
				//发送稀疏数据
				MPI_Isend(value_buffer,nnz * (isize + vsize), MPI_CHAR,right, MessageType::kScatterReduce,comm,&requests[0]);
			}else{
				MPI_Isend(base, block_size * vsize, MPI_CHAR, right, MessageType::kScatterReduce, comm,&requests[0]);
			}
		}else{
			MPI_Isend(base, block_size * vsize, MPI_CHAR, right, MessageType::kScatterReduce, comm,&requests[0]);
		}
		 base=buffer+blocks[my_index].first;
		 block_size = blocks[my_index].second;
 
		MPI_Irecv(recv_buffer, block_size * vsize, MPI_CHAR, left, MessageType::kScatterReduce, comm, &requests[1]);
        MPI_Waitall(2, requests, statuses);
		int nnz=0;
		MPI_Get_count(&statuses[1], MPI_CHAR, &nnz);
		 
		if(nnz<block_size * vsize)//接受到稀疏数据
		{
		    CHECK_EQ(nnz % (isize + vsize), 0);
			nnz = nnz / (isize + vsize);
			memcpy(index_buffer, recv_buffer + nnz, nnz * isize);
			for(int j=0;j<nnz;j++)
			{
				//base[index_buffer[j]]+=recv_buffer[j];
				Reduce<O>(base[index_buffer[j]], recv_buffer[j]);
			}
		}else{//接收到稠密数据
		   for (int j = 0; j < block_size; ++j) {
               // base[j]+=recv_buffer[j];
				 Reduce<O>(base[j], recv_buffer[j]);
            }
		}
	 
	}
  
	//algather阶段
	need_check=true;
	int total_num=1;
	int number=pow(2,total_num);
	while(number<worker_number)
	{
		total_num++;
		number=pow(2,total_num);
	}
	
	//cout<<"********"<<endl;
	int sendBlockSize= blocks[my_index].second;
	T *base=buffer+blocks[my_index].first;
	int start_index=blocks[my_index].first;
	int recvBlockSize=0;
	int dis=1; 
	int thred=2;
	T *sBuf=new T[count];
	int *index=new int[count];
	int nnz_cnt=0;
	for(int i=0;i<sendBlockSize;i++)
	{
		if(base[i]!=0)
		{
			sBuf[nnz_cnt]=base[i];
			index[nnz_cnt]=i+start_index;
			nnz_cnt++;
		}
		//sBuf[i]=base[i];
	}
	 
	T *send_buf=new T[count*2];
	T *recv_buf=new T[count];
	//for(int i=0;i<total_num;i++)
	for(int i=total_num;i>0;i--)
	{
		MPI_Status statuses[2];
        MPI_Request requests[2];
		//dis=pow(2,i);
		dis=pow(2,i-1);
	    thred=dis*2;
		 
	    if(my_index%thred+dis<thred)
	    {
		    right=worker_list[my_index+dis];
	    }else{
		
		    right=worker_list[my_index-dis];
	    }
        
		memcpy(send_buf,sBuf,nnz_cnt*vsize); 
		memcpy(send_buf+nnz_cnt,index,nnz_cnt*isize);
		
		MPI_Isend(send_buf,nnz_cnt*(vsize+isize),MPI_CHAR,right,MessageType::kAllGather,comm,&requests[0]);
		
		int rcv_acnt=0;
		MPI_Status recvstatus;
		MPI_Probe(right,MPI_ANY_TAG,comm,&recvstatus);
		MPI_Get_count(&recvstatus, MPI_CHAR, &rcv_acnt);
		char *recv_buf=new char[rcv_acnt];
		MPI_Irecv(recv_buf,rcv_acnt,MPI_CHAR,right,MessageType::kAllGather,comm,&requests[1]);
		MPI_Wait(&requests[0],&statuses[0]); 
        MPI_Wait(&requests[1],&statuses[1]);
		
		recvBlockSize=rcv_acnt/(vsize+isize);
		memcpy(sBuf+nnz_cnt,recv_buf,recvBlockSize*vsize);
		memcpy(index+nnz_cnt,recv_buf+recvBlockSize*vsize,recvBlockSize*isize);
		nnz_cnt+=recvBlockSize;
		delete[] recv_buf;
	}
 	 
	for(int i=0;i<nnz_cnt;i++)
	{
		 
		  buffer[index[i]]=sBuf[i];
	}
	 
	delete[] send_buf;
	delete[] index;
	delete[] sBuf;
	delete[] recv_buffer;
    delete[] send_buffer;
    delete[] index_buffer;
}

//split-allgather的算法实现
template<class O, class T>
void SplitAlgather(T *buffer, int count, int id, std::vector<int> &worker_list, MPI_Comm comm) {
int worker_number = worker_list.size();
    if (worker_number == 1) return;
    if (count < 512 * worker_number) {
        SimpleAllreduce<O>(buffer, count, id, worker_list, comm);
        return;
    }

    
	int block_size=count/worker_number;
	std::map<int,std::pair<int,int>> blocks; //记录每块的起始地址和块大小
	for(int i=0;i<worker_number;i++)
    {
		blocks[i].first=i*block_size;
	    blocks[i].second=block_size;
	}		
	blocks[worker_number-1].second = count - (worker_number - 1) * block_size;
	int right=-1;
	int left=-1;
	int my_index=-1;
	
	for(int i=0;i<worker_number;i++)
	{
		if(worker_list[i]==id)
		{
			my_index=i;
			break;
		}
	}
	 
	bool need_check = true;
	int isize = sizeof(int);
    int vsize = sizeof(T);
	T *send_buffer=new T[count];
	T *value_buffer = new T[count];
	int *index_buffer=new int[count];
	 
	T *recv_buffer=new T[count];
	//MPI_Request req[2];
	//MPI_Status status[2];
 
	//split 阶段
    for(int i=0;i<worker_number-1;++i)
	{
		left=worker_list[(my_index-i-1+worker_number)%worker_number];
	 	right=worker_list[(my_index+i+1)%worker_number];
		MPI_Status statuses[2];
        MPI_Request requests[2];
		T *base=buffer+blocks[right].first;
		block_size=blocks[right].second;
 
		if(need_check)//是否过滤标志
		{
			int nnz=0;
			for(int j=0;j<block_size;j++)
			{
				if(base[j]!=0)
				{
					++nnz;
				}
			}
			 
			if(nnz*(isize+vsize)<block_size * vsize)//满足稀疏传送条件
			{
				int k=0;
				for(int j=0;j<block_size;++j)
				{
					if(base[j]!=0)
					{
						value_buffer[k]=base[j];
						index_buffer[k]=j;
						++k;
					}					
				}
				memcpy(value_buffer+nnz,index_buffer,isize*nnz); //将对应索引值放在value值的后面
				//发送稀疏数据
				MPI_Isend(value_buffer,nnz * (isize + vsize), MPI_CHAR,right, MessageType::kScatterReduce,comm,&requests[0]);
			}else{
				MPI_Isend(base, block_size * vsize, MPI_CHAR, right, MessageType::kScatterReduce, comm,&requests[0]);
			}
		}else{
			MPI_Isend(base, block_size * vsize, MPI_CHAR, right, MessageType::kScatterReduce, comm,&requests[0]);
		}
		 base=buffer+blocks[my_index].first;
		 block_size = blocks[my_index].second;
 
		MPI_Irecv(recv_buffer, block_size * vsize, MPI_CHAR, left, MessageType::kScatterReduce, comm, &requests[1]);
        MPI_Waitall(2, requests, statuses);
		int nnz=0;
		MPI_Get_count(&statuses[1], MPI_CHAR, &nnz);
		 
		if(nnz<block_size * vsize)//接受到稀疏数据
		{
		    CHECK_EQ(nnz % (isize + vsize), 0);
			nnz = nnz / (isize + vsize);
			memcpy(index_buffer, recv_buffer + nnz, nnz * isize);
			for(int j=0;j<nnz;j++)
			{
				//base[index_buffer[j]]+=recv_buffer[j];
				Reduce<O>(base[index_buffer[j]], recv_buffer[j]);
			}
		}else{//接收到稠密数据
		   for (int j = 0; j < block_size; ++j) {
               // base[j]+=recv_buffer[j];
				 Reduce<O>(base[j], recv_buffer[j]);
            }
		}
	 
	}
  
	//algather阶段
	need_check=true;
	int total_num=1;
	int number=pow(2,total_num);
	while(number<worker_number)
	{
		total_num++;
		number=pow(2,total_num);
	}
	
	//cout<<"********"<<endl;
	int sendBlockSize= blocks[my_index].second;
	T *base=buffer+blocks[my_index].first;
	int start_index=blocks[my_index].first;
	int recvBlockSize=0;
	int dis=1; 
	int thred=2;
	T *sBuf=new T[count];
	int *index=new int[count];
	int nnz_cnt=0;
	for(int i=0;i<sendBlockSize;i++)
	{
		if(base[i]!=0)
		{
			sBuf[nnz_cnt]=base[i];
			index[nnz_cnt]=i+start_index;
			nnz_cnt++;
		}
		//sBuf[i]=base[i];
	}
	 
	T *send_buf=new T[count*2];
	T *recv_buf=new T[count];
	for(int i=0;i<total_num;i++)
	{
		MPI_Status statuses[2];
        MPI_Request requests[2];
		dis=pow(2,i);
	    thred=dis*2;
		 
	    if(my_index%thred+dis<thred)
	    {
		    right=worker_list[my_index+dis];
	    }else{
		
		    right=worker_list[my_index-dis];
	    }
        
		memcpy(send_buf,sBuf,nnz_cnt*vsize); 
		memcpy(send_buf+nnz_cnt,index,nnz_cnt*isize);
		
		MPI_Isend(send_buf,nnz_cnt*(vsize+isize),MPI_CHAR,right,MessageType::kAllGather,comm,&requests[0]);
		
		int rcv_acnt=0;
		MPI_Status recvstatus;
		MPI_Probe(right,MPI_ANY_TAG,comm,&recvstatus);
		MPI_Get_count(&recvstatus, MPI_CHAR, &rcv_acnt);
		char *recv_buf=new char[rcv_acnt];
		MPI_Irecv(recv_buf,rcv_acnt,MPI_CHAR,right,MessageType::kAllGather,comm,&requests[1]);
		MPI_Wait(&requests[0],&statuses[0]); 
        MPI_Wait(&requests[1],&statuses[1]);
		
		recvBlockSize=rcv_acnt/(vsize+isize);
		memcpy(sBuf+nnz_cnt,recv_buf,recvBlockSize*vsize);
		memcpy(index+nnz_cnt,recv_buf+recvBlockSize*vsize,recvBlockSize*isize);
		nnz_cnt+=recvBlockSize;
		delete[] recv_buf;
	}
 	 
	for(int i=0;i<nnz_cnt;i++)
	{
		 
		  buffer[index[i]]=sBuf[i];
	}
	 
	delete[] send_buf;
	delete[] index;
	delete[] sBuf;
	delete[] recv_buffer;
    delete[] send_buffer;
    delete[] index_buffer;
}

 //split-scannter的算法实现
template<class O, class T>
void SplitScannterAllreduce(T *buffer, int count, int id, std::vector<int> &worker_list, MPI_Comm comm) {
    int worker_number = worker_list.size();
    if (worker_number == 1) return;
    if (count < 512 * worker_number) {
        SimpleAllreduce<O>(buffer, count, id, worker_list, comm);
        return;
    }

    int block_size = count / worker_number;
    std::map<int, std::pair<int, int>> blocks;
    for (int i = 0; i < worker_number; ++i) {
        blocks[i].first = i * block_size;
        blocks[i].second = block_size;
    }
    blocks[worker_number - 1].second = (count - (worker_number - 1) * block_size);
    int left = -1, right = -1, my_index = -1;
    for (int i = 0; i < worker_number; ++i) {
        if (worker_list[i] == id) {
            my_index = i;
          //  left = worker_list[(i - 1 + worker_number) % worker_number];
           // right = worker_list[(i + 1) % worker_number];
            break;
        }
    }
   // CHECK(left != -1 && right != -1 && my_index != -1);

    bool need_check = true;
    int isize = sizeof(int);
    int vsize = sizeof(T);
    int *index_buffer = new int[count];
    T *value_buffer = new T[count];
    T *recv_buffer = new T[count];
   
    for (int i = 0; i < worker_number - 1; ++i) {
		left=worker_list[(my_index-i-1+worker_number)%worker_number];
		right=worker_list[(my_index+i+1)%worker_number];
        MPI_Status statuses[2];
        MPI_Request requests[2];
        //T *base = buffer + blocks[right].first;
        block_size = blocks[right].second;
		T *base=new T[block_size];
		int start=blocks[right].first;
		int nnz = 0;
        for (int j = 0; j < block_size; ++j)
		{
			base[j]=buffer[j+start];
            if (base[j] != 0) 
			{
                ++nnz;
            }
        }
        //如果上一次迭代接收到了稠密的块，那么本次传输就不用检查了，必然稠密
        if (need_check) {
           /* int nnz = 0;
            for (int j = 0; j < block_size; ++j) {
                if (base[j] != 0) {
                    ++nnz;
                }
            }*/
            //如果满足稀疏传输条件则采用稀疏传输
            if (nnz * (isize + vsize) < block_size * vsize) {
                int k = 0;
                for (int j = 0; j < block_size; ++j) {
                    if (base[j] != 0) {
                        value_buffer[k] = base[j];
                        index_buffer[k] = j;
                        ++k;
                    }
                }
                //把元素的值和元素的索引放到同一个内存空间上，一同发送，不分两次发送了
                memcpy(value_buffer + nnz, index_buffer, isize * nnz);
                MPI_Isend(value_buffer, nnz * (isize + vsize), MPI_CHAR, right,
                          MessageType::kScatterReduce, comm, &requests[0]);
            } else {
                MPI_Isend(base, block_size * vsize, MPI_CHAR, right, MessageType::kScatterReduce, comm,
                          &requests[0]);
            }
        } else {
            MPI_Isend(base, block_size * vsize, MPI_CHAR, right, MessageType::kScatterReduce, comm,
                      &requests[0]);
        }
        //base = buffer + blocks[my_index].first;
        block_size = blocks[my_index].second;
		start= blocks[my_index].first;
        MPI_Irecv(recv_buffer, block_size * vsize, MPI_CHAR, left, MessageType::kScatterReduce, comm, &requests[1]);
        MPI_Waitall(2, requests, statuses);
        nnz = 0;
        MPI_Get_count(&statuses[1], MPI_CHAR, &nnz);
        //如果接收到的字节数少于这个块应有的大小，说明左边进程发送的是稀疏数据
        if (nnz < block_size * vsize) {
            CHECK_EQ(nnz % (isize + vsize), 0);
            nnz = nnz / (isize + vsize);
            //分离值和索引
            memcpy(index_buffer, recv_buffer + nnz, nnz * isize);
            for (int j = 0; j < nnz; ++j) {
               // Reduce<O>(base[index_buffer[j]], recv_buffer[j]);
				buffer[index_buffer[j]+start]+=recv_buffer[j];
            }
           // need_check = true;
        } else {
            CHECK_EQ(nnz, block_size * vsize);
            for (int j = 0; j < block_size; ++j) {
               // Reduce<O>(base[j], recv_buffer[j]);
				buffer[j+start]+=recv_buffer[j];
            }
           // need_check = false;
        }
         delete []base;
        
    }
    need_check = true;
    for (int i = 0; i < worker_number - 1; ++i) {
        MPI_Status statuses[2];
        MPI_Request requests[2];
        
		right=worker_list[(my_index+1)%worker_number];
		left=worker_list[(my_index+worker_number-1)%worker_number];
		int sendBlockId=worker_list[(my_index+worker_number-i)%worker_number];
		int recvBlockId=worker_list[(my_index+worker_number-i-1)%worker_number];
		//T *base = buffer + blocks[sendBlockId].first;
        block_size = blocks[sendBlockId].second;
		int start=blocks[sendBlockId].first;
		T *base=new T[block_size];
		int nnz = 0;
        for (int j = 0; j < block_size; ++j) {
			base[j]=buffer[start+j];
            if (base[j] != 0) {
                    ++nnz;
            }
        }
        if (need_check) {
            
            if (nnz * (isize + vsize) < block_size * vsize) {
                int k = 0;
                for (int j = 0; j < block_size; ++j) {
                    if (base[j] != 0) {
                        value_buffer[k] = base[j];
                        index_buffer[k] = j;
                        ++k;
                    }
                }
                memcpy(value_buffer + nnz, index_buffer, isize * nnz);
                MPI_Isend(value_buffer, nnz * (isize + vsize), MPI_CHAR, right,
                          MessageType::kAllGather, comm, &requests[0]);
            } else {
                MPI_Isend(base, block_size * vsize, MPI_CHAR, right, MessageType::kAllGather, comm,
                          &requests[0]);
            }
        } else {
            MPI_Isend(base, block_size * vsize, MPI_CHAR, right, MessageType::kAllGather, comm,
                      &requests[0]);
        }
        //base = buffer + blocks[recvBlockId].first;
        block_size = blocks[recvBlockId].second;
		start = blocks[recvBlockId].first;
        MPI_Irecv(recv_buffer, block_size * vsize, MPI_CHAR, left, MessageType::kAllGather, comm, &requests[1]);
        MPI_Waitall(2, requests, statuses);
        nnz = 0;
        MPI_Get_count(&statuses[1], MPI_CHAR, &nnz);
        if (nnz < block_size * vsize) {
            CHECK_EQ(nnz % (isize + vsize), 0);
            nnz = nnz / (isize + vsize);
            memcpy(index_buffer, recv_buffer + nnz, nnz * isize);
          
            for (int j = 0; j < nnz; ++j) {
                //base[index_buffer[j]] = recv_buffer[j];
				buffer[index_buffer[j]+start] = recv_buffer[j];
            }
           // need_check = true;
        } else {
            CHECK_EQ(nnz, block_size * vsize);
            for (int j = 0; j < block_size; ++j) {
                //base[j] = recv_buffer[j];
				buffer[j+start] = recv_buffer[j];
            }
            //need_check = false;
        }
         delete []base;
        //send_block_index = recv_block_index;
       // recv_block_index = (recv_block_index - 1 + worker_number) % worker_number;
    }

    delete[] recv_buffer;
    delete[] value_buffer;
    delete[] index_buffer;
}

}

#endif //SPARSEALLREDUCE_RING_ALLREDUCE_H
