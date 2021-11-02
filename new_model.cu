//////////////////////////////////////////New code////////////////////////////////////

#include<iostream>
#include<vector>
#include<set>
#include<fstream>
#include <chrono>
#include<stdio.h>


using namespace std::chrono;
using std::ifstream;
ifstream indata ; 

#define ll long long int 


using namespace std;
/*----------------------------------------------------------------------------*/

__device__ bool check_edge(ll* edges  , ll left ,ll right , ll ff , ll ss )
{
 
    while(left<=right)
    {
        ll center= (right-left)/2 + left ; 
   //  printf("%d %d %d %d \n" , edges[center*2] , edges[center*2+1] , ff , ss  ) ; 
        if(edges[center*2]== ff && edges[center*2+1] == ss)
        return 1 ;  //Found 

        if(edges[center*2] > ff )
        {
            right = center -1 ; 
        }else if(edges[center*2] < ff )
        left = center+1 ; 
        else if (edges[center*2+1] < ss)
        left = center +1 ; 
        else 
        right = center -1 ; 
     }
    return 0 ; 
}

/*----------------------------------------------------------------------------*/

__device__ int match_value(ll* last_cliq ,ll center ,ll index ,ll cli )
{
    for(ll i =0 ;i <cli-1 ;i++)
    {
        if(last_cliq[center*cli+i] > last_cliq[index*cli+i+1])
        return 1 ;
        else if (last_cliq[center*cli+i] < last_cliq[index*cli+i+1])
        return -1 ; 
    }
    return 0 ; 
}
/*----------------------------------------------------------------------------*/
__device__ ll find_starting_point(ll* last_cliq ,ll left ,  ll right  ,ll cli , ll index )
{
    ll ans = -1 ; 
    while(left <= right )
    {
        ll center = (right-left)/2 +left ;
        if(match_value(last_cliq , center , index , cli ) == 0 )
        {
            ans = center ; 
            right =center -1 ; 
        } else if (match_value(last_cliq , center , index , cli ) == -1 )
        {   
            
            
            //center is less than index 
            left = center +1 ; 
        }else 
          right = center -1 ; 
    }
 return  ans ; 
} 
/*----------------------------------------------------------------------------*/
__global__ void count_next_cli (ll* edges ,ll* last_cliq ,ll n ,ll m ,ll current_cli_size ,ll*  start_3_cli , ll cli  )
{
    
      ll index = blockDim.x * blockIdx.x + threadIdx.x ;
      if(index < current_cli_size)
      {   
           
           start_3_cli[index] = 0 ; //Preset value 
          ll starting_point = find_starting_point(last_cliq ,0 ,  current_cli_size-1 , cli , index ) ;
          if(starting_point == -1 )
          return ;    
          for( ; starting_point < current_cli_size && match_value(last_cliq , starting_point , index , cli )==0; starting_point++)
          {   bool added = 1 ; 

              //printf("%d\n" , starting_point ) ; 
              //So we need to see that is the last element connected to all 
              for(ll j = 0 ; j <cli ; j++ )
              {
                  if(check_edge(edges , 0 , m -1 , last_cliq[starting_point*cli+cli-1] , last_cliq[index*cli+j]  )==0)
                  {
                      added = 0 ;
                  }
              }
           start_3_cli[index] += added  ; 
                  
          }
      }
}


__global__ void make_next_cli (ll* edges ,ll* last_cliq ,ll n ,ll m ,ll current_cli_size ,ll*  places , ll cli , ll* update  )
{
      ll index = blockDim.x * blockIdx.x + threadIdx.x ;
      if(index < current_cli_size)
      {   
          ll starting_point = find_starting_point(last_cliq ,0 ,  current_cli_size-1 , cli , index ) ;
          if(starting_point == -1 )
          return ; 
          ll writing_ptr = places[index] ; 
          for( ; starting_point < current_cli_size && match_value(last_cliq , starting_point , index , cli )==0; starting_point++)
          {   bool added = 1 ; 
            
              //So we need to see that is the last element connected to all 
              for(ll j = 0 ; j <cli ; j++ )
              {
                  if(check_edge(edges , 0 , m -1 , last_cliq[starting_point*cli+cli-1] , last_cliq[index*cli+j]  )==0)
                  added = 0 ;

                    
              }
           if(added)
            {
                for(ll i =0 ;i < cli  ; i++)
                {
                    update[writing_ptr*(cli+1)+i] = last_cliq[index*cli+i] ; 
                 
                }
                    update[writing_ptr*(cli+1)+cli] = last_cliq[starting_point*cli+cli-1] ; //I want the last element 
                  writing_ptr++ ; 
            }
                  

           
          }
      }
}





void  make_two_cli(ll m  , vector<ll>& host_edges , vector<ll>& host_two_cli  )
{ 
    set<pair<ll,ll>> tem ;  
    for(ll i =0;i<m;i++)
    {
        ll a,b ; indata>>a>>b ; 
        tem.insert({a,b})  ;     
        tem.insert({b,a})  ;    //Considering the undirected edges 
    
    }
    for(auto i : tem)
    {
        /* Directed 
            if(i.first < i.second && (tem.find({i.second , i.first }) != tem.end()) )
                {host_two_cli.push_back(i.first) ;  host_two_cli.push_back(i.second) ;} 
            host_edges.push_back(i.first);host_edges.push_back(i.second); //It will be automatically sorted 
     
        */
        
        /*Undirected*/
        
        host_edges.push_back(i.first);host_edges.push_back(i.second); //It will be automatically sorted 
        
      if(i.first < i.second )
         {host_two_cli.push_back(i.first) ;  host_two_cli.push_back(i.second) ; }
       
     
    }
 //So I have updated  the value of both 
}

ll find_kcliq(ll k)
{
    ll  m ;indata>>m ;
    ll n = 1 ;
    /// take the input 
    vector<ll> host_edges; 
    vector<ll> host_cli;
    make_two_cli(m , host_edges , host_cli) ;
    ll current_cli_size = host_cli.size()/2 ; //Num of cliques  
    for(auto i : host_cli)
      n = max(i,n) ; 
    if(k==1)
    return n ; 
    if(k==2)
    return current_cli_size ; 
    n++ ; 
    m*=2 ; //only for undirected graph 
    ll cli = 2 ;
 
    //Upload the edges data on the 

       ll* device_edges   ; 
       cudaMalloc( &device_edges , sizeof(ll)*m*2) ;
       cudaMemcpy(device_edges ,host_edges.data()  , sizeof(ll)*(m*2) , cudaMemcpyHostToDevice )  ;   
       ll* last_cliq ;
       cudaMalloc( &last_cliq , sizeof(ll)*current_cli_size*cli) ;
       cudaMemcpy(last_cliq, host_cli.data()  , sizeof(ll)*current_cli_size*cli , cudaMemcpyHostToDevice )  ;
      //  cout << current_cli_size <<"\n" ; 
        while(cli<k)
        {
            //Find the number of cliques for cli+1 
            ll* device_number_cli ; 
            cudaMalloc( &device_number_cli , sizeof(ll)*current_cli_size) ;
            int threadsPerBlock = 100 ;
            ll blocksPerGrid = ( (current_cli_size)+ threadsPerBlock - 1) / threadsPerBlock;
            count_next_cli<<<blocksPerGrid, threadsPerBlock>>> (device_edges , last_cliq ,n , m , current_cli_size ,  device_number_cli , cli );
             
            ll new_next_cliq = 0 ; 
            vector<ll> number_cli(current_cli_size);
            cudaMemcpy(number_cli.data() , device_number_cli, sizeof(ll)*(current_cli_size) , cudaMemcpyDeviceToHost);
            for(ll i=0;i<current_cli_size ; i++)
            {
                ll tem = number_cli[i] ; 
                number_cli[i] = new_next_cliq ; //It becomes the starting point  
                new_next_cliq +=  tem ; 
             
            }
         
            if(cli+1 < k)
            {
                //So I would require to add new value 
                ll*  dummy_last_cli ; 
                cudaMalloc(&dummy_last_cli , sizeof(ll)*new_next_cliq*(cli+1)) ; //So this will be the length 
                //We would had to copy it back 
                cudaMemcpy(  device_number_cli,number_cli.data() ,  sizeof(ll)*(current_cli_size) , cudaMemcpyHostToDevice);
                make_next_cli<<<blocksPerGrid, threadsPerBlock>>> (device_edges , last_cliq ,n , m , current_cli_size ,  device_number_cli , cli , dummy_last_cli );

                cudaFree(last_cliq) ; 
                last_cliq =  dummy_last_cli ;        
                vector<ll> ads (new_next_cliq*(cli+1 )) ; 

                cudaMemcpy(ads.data() ,last_cliq, sizeof(ll)*(new_next_cliq*(cli+1 )) , cudaMemcpyDeviceToHost);


                //We need no memory copy 
            }

            //We need to update the values now 
            cli++ ; 
            current_cli_size = new_next_cliq ; 



        }   
 
      cudaFree(last_cliq) ; 
      cudaFree(device_edges) ; 
      return current_cli_size ; 
       
 
}
int main()
{
    
    ll k  ;
    string file_name   ;
  
    cin>>file_name >> k ; 
//    cout<<file_name <<" " << k << "\n" ; 
    indata.open(file_name); // opens the file
    if(!indata) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
    }
  milliseconds ms = duration_cast< milliseconds >(
    system_clock::now().time_since_epoch()
);
 
    ll ms1 =std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ll ans  =  find_kcliq(k) ;
    ll ms2 = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    
    cout<< ans  <<"\nExecution Time: "<< ms2 - ms1 <<"ms\n" ;

}
