//////////////////////////////////////////New code////////////////////////////////////

#include<iostream>
#include<vector>
#include<set>
#include<fstream>
#include <chrono>


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

__global__ void make_count (ll* device_edges , ll* device_last_ele  ,ll current_cli_size ,ll m ,ll* device_count )
{
    
      ll index = blockDim.x * blockIdx.x + threadIdx.x ;
    if(index < current_cli_size )
    {
        device_count[index] =0  ; 
        for(ll i = index +1 ; i<current_cli_size && device_last_ele[index*2] == device_last_ele[i*2] ; i++)
        {
            if(check_edge(device_edges  , 0 , m-1 , device_last_ele[index*2+1] , device_last_ele[i*2+1]))
            device_count[index]++ ; //Since it is sorted we need not worry about anything other  
        }
    }
}
    
/*----------------------------------------------------------------------------*/

__global__ void  make_next(ll* device_edges ,ll* device_last_ele  ,ll current_cli_size ,ll m ,ll* device_count ,ll* device_new_last )
{
    
      ll index = blockDim.x * blockIdx.x + threadIdx.x ;
    if(index < current_cli_size )
    {
        ll starting_ele = device_count[index] ;
        for(ll i = index+ 1  ; i<current_cli_size && device_last_ele[index*2] == device_last_ele[i*2] ; i++)
        {
            if(check_edge(device_edges  , 0 , m-1 , device_last_ele[index*2+1] , device_last_ele[i*2+1]))
            {
                device_new_last[2*starting_ele] = index ; 
                device_new_last[1+2*starting_ele] = device_last_ele[i*2+1] ; 
                starting_ele++ ; 
            }
            
        }
    }
}
       
/*----------------------------------------------------------------------------*/







void  make_two_cli(ll m  , vector<ll>& host_edges , vector<ll>& host_two_cli  )
{ 
    set<pair<ll,ll>> tem ;  
    for(ll i =0;i<m;i++)
    {
        ll a,b ; indata>>a>>b ; 
        tem.insert({a,b})  ;    //Considering the directed edges 
        tem.insert({b,a})  ;    //Considering the directed edges 
    
    }
    for(auto i : tem)
    {
        
        /*
        //Directed
      if(i.first < i.second && (tem.find({i.second , i.first }) != tem.end()) )
         {host_two_cli.push_back(i.first) ;  host_two_cli.push_back(i.second) ;} 
        */
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
    
    m =  host_edges.size()/2  ; 
    
    ll cli = 2 ;
    ll *device_last_ele  , *device_edges ; 
    
    
    cudaMalloc( &device_edges , sizeof(ll)*m*2) ;
    cudaMalloc( &device_last_ele , sizeof(ll)*current_cli_size*2) ;

    cudaMemcpy(device_edges ,host_edges.data()  , sizeof(ll)*(m*2) , cudaMemcpyHostToDevice )  ;   
    cudaMemcpy(device_last_ele ,host_cli.data()  , sizeof(ll)*(current_cli_size*2) , cudaMemcpyHostToDevice )  ;   
    
    while(cli < k )
    {
        //We need to count the number of next nodes 

        ll *device_count ; 
        vector<ll> count(current_cli_size) ; 
        cudaMalloc( &device_count  , sizeof(ll)*current_cli_size) ;
        

        int threadsPerBlock = 100 ;
        ll blocksPerGrid = ( (current_cli_size)+ threadsPerBlock - 1) / threadsPerBlock;
        make_count<<<blocksPerGrid, threadsPerBlock>>> (device_edges , device_last_ele  , current_cli_size , m , device_count );
        cudaMemcpy(count.data() , device_count, sizeof(ll)*(current_cli_size) , cudaMemcpyDeviceToHost);

        ll new_cli_size =0  ;
        for(ll i =0;i<current_cli_size ; i++)
        {
            ll temp =  count[i] ; 
            count[i] = new_cli_size ; 
         new_cli_size += temp ; 
        }
      
        cudaMemcpy(device_count,count.data()   , sizeof(ll)*current_cli_size, cudaMemcpyHostToDevice )  ;
  
     if(new_cli_size==0)
     {
         current_cli_size = 0 ; 
      break ; 
     }
     
     if(1+cli < k   )
     {
         //Updete the value
        ll  *device_new_last ; 
        cudaMalloc( &device_new_last  , sizeof(ll)*new_cli_size*2) ;

        make_next<<<blocksPerGrid, threadsPerBlock>>> (device_edges , device_last_ele  , current_cli_size , m , device_count  , device_new_last );
        cudaFree(device_last_ele) ; 
        device_last_ele = device_new_last ; 
       
     }
     

     cli++ ; 
     current_cli_size = new_cli_size ; 
 
    }
 
    cudaFree(device_edges) ; 
    cudaFree(device_last_ele) ;
    return  current_cli_size ; 
 
    


}
int main()
{
    
    ll k ;
    //cin>> k ;
    string file_name   ;
    cin>>file_name >>k ; 
    
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
    
    cout<< ans  <<"\nExecution Time: "<< ms2 - ms1 <<"ms" ;

}
