#include<iostream>
#include<vector>
#include<algorithm>
#include <ctime>
#include<map>

using namespace std;
#define ll long long int 



#include <fstream>
using std::ifstream;

     ifstream indata; // indata is like cin



//The vector file required 

template<typename T>
class LocalVector
{
private:
    T* m_begin;
    T* m_end;

    size_t capacity;
    __device__ void expand() {
        capacity *= 2;
        size_t tempLength = (m_end - m_begin);
        T* tempBegin = new T[capacity];

        memcpy(tempBegin, m_begin, tempLength * sizeof(T));
        delete[] m_begin;
        m_begin = tempBegin;
        m_end = m_begin + tempLength;
        length = static_cast<size_t>(m_end - m_begin);
    }
public:
    size_t length;
 
    __device__  explicit LocalVector() : length(0), capacity(16) {
        m_begin = new T[capacity];
        m_end = m_begin;
    }
    __device__ T& operator[] (unsigned int index) {
        return *(m_begin + index);//*(begin+index)
    }
    __device__ T* begin() {
        return m_begin;
    }
    __device__ T* end() {
        return m_end;
    }
    __device__ ~LocalVector()
    {
        delete[] m_begin;
        m_begin = nullptr;
    }

    __device__ void add(T t) {

        if ((m_end - m_begin) >= capacity) {
            expand();
        }

        new (m_end) T(t);
        m_end++;
        length++;
    }
    __device__ T pop() {
        T endElement = (*m_end);
        delete m_end;
        m_end--;
        return endElement;
    }

    __device__ size_t getSize() {
        return length;
    }
};



/////////////////////////////////



//Sorting algo
void sort( ll *edges , ll m )
{
  vector<pair<ll,ll>> arr ;
  for(ll i =0;i<m;i++)
  arr.push_back({edges[i*2] , edges[i*2+1]}) ; 
  sort(arr.begin() , arr.end()) ; 
 for(ll i =0;i<m;i++)
 {
     edges[i*2] = arr[i].first ; 
     edges[i*2+1] = arr[i].second ; 
 }

}

///////////////////////////////

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

__device__ ll first_node_greater_than_n(ll* edges , ll left ,ll right , ll n )
{   ll ans = -1  ; 
    while(left<= right )
    {
        ll center = (right - left)/2 + left ; 

        if(edges[center*2] > n )
        {
            ans =  center ; 
            right = center-1; 
        }else 
          left = center +1 ; 
    }
 return ans ; 

}


//use two clique beacause its distinct 
__device__ ll start_ele(ll * two_cliq , ll m , ll ele )
{
    ll left = 0 , right = m-1 ; 
    ll ans = -1 ; 

 while(left <= right )
 {
    ll center = (right - left ) /2 + left ; 
  if(two_cliq[center*2] == ele  ) 
  ans = center ; 

  if(two_cliq[center*2] >=  ele )
  right = center -1 ; 
  else 
  left = center +1 ;  
 
 }
 return ans ; 

}




//---------------------------------------------------------------------------------------------------------
//For k<=4 

__device__ ll sol(ll*edges , LocalVector<ll> &arr  , ll k , ll now ,  ll n ,ll m , ll* two_cliques )
{
 
    if(k==now)
    return 1 ; 
    ll ans =0 ; 

    if(k-2 >= now )
    {   
        ll starting_index = arr[now-1] ; 
        starting_index  = first_node_greater_than_n(two_cliques , 0, (m-1)/2,  starting_index) ; 
        if(starting_index == -1 )
        return 0 ; 
  //   printf("%d\n"  , starting_index ) ;

        for( ; starting_index < m/2 ;starting_index++)
        {
            bool to_add = 1 ;
         for(ll i = 0 ; i < now ;i++)
         {
             if( !( check_edge(edges , 0 , m-1 , arr[i] , two_cliques[starting_index*2]   ) && check_edge(edges , 0 , m-1 , arr[i] , two_cliques[starting_index*2+1]   )  || arr[i]  >= two_cliques[starting_index*2]) )
            {
                to_add = 0 ;
                break ; 
            }
         }
          if(to_add == 1 )
          {   
              arr[now] = two_cliques[starting_index*2] ;

              arr[now+1] = two_cliques[starting_index*2+1] ; 
              
              ans +=  sol(edges , arr , k , now+2 , n , m , two_cliques) ; 
              
          }
        }
       
    }else 
    {
   
    ll starting_element =  start_ele(two_cliques , m/2 , arr[now-1]  ) ; 
  //  printf("%d\n", starting_element ) ; 
      
        if(starting_element == -1 )
        return 0 ; 
     ll sd = arr[3] - arr[2] ; 
   //     printf("%d\n" , sd ) ;
      for(    ; starting_element < m/2 && two_cliques[starting_element*2] == arr[now-1] ; starting_element++  )
      {
           ll i = two_cliques[starting_element*2+1] ;
    
           bool added = 1 ; 
          for(ll j = 0 ; j < now ; j++ )
          {
              if(check_edge(edges , 0 , m-1,  i , arr[j]) == 0 )
              {
                  added = 0 ;
                  break ;  
              }
          }
       //So we can work with i 
        if(added)
        {
       //     printf("%d %d %d %d %d %d %d %d %d %d\n" , arr[0] , arr[0] , arr[1] , arr[1] , arr[2] , arr[2] ,arr[3] ,arr[3] , arr[4] ,arr[4]) ;
            arr[now] = i ; 
            ans +=  1 ; 
        }
       } 
    }
 return ans ; 
}

__global__ void doi( ll* edges , ll * updates , ll *two_cliques , ll n , ll m , ll k )
{
    
    ll index = blockDim.x * blockIdx.x + threadIdx.x ;
    if(index <  m/2)
    {
        
     
        if(k==2)
        {
            updates[index] = 1 ; 
            return ; 
        }
      
         //Now we will find the 3-cliques 

        ll current_size = 0 ;

     
              //  printf("%d\n", two_cliques[2*index+1]+1)  ;
                   LocalVector<ll> dum  ; 
                    dum.add( two_cliques[index*2]  ) ;
                    dum.add( two_cliques[index*2+1]  ) ;
                  for(ll i =0;i<k-2 ; i++ )
                  dum.add(0) ;
            current_size += sol(edges , dum , k , 2 , n , m , two_cliques ) ; 
         

     
         updates[index]  = current_size  ; 
         return ; 
    
    }    
}




//--------------------------------------------------------------------------------

__global__ void count_three_cliques (ll* edges ,ll*  two_cliques,ll    n ,ll  m , ll* start_3_cli    ) 
      {
          ll index = blockDim.x * blockIdx.x + threadIdx.x ;
       if(index<m/2)
       {
       
          ll ans = 0 ;    
          ll starting_element =  start_ele(two_cliques , m/2 , two_cliques[index*2+1]  ) ; 

          if(starting_element == -1 )
          return  ; 
          for(    ; starting_element < m/2 && two_cliques[starting_element*2] ==  two_cliques[index*2+1] ; starting_element++  )
          {
           ll i = two_cliques[starting_element*2+1] ;
    
           bool added = 1 ; 
          for(ll j = 0 ; j < 2; j++ )
          {
              if(check_edge(edges , 0 , m-1,  i , two_cliques[index*2+j] ) == 0 )
              {
                  added = 0 ;
                  break ;  
              }
          }
       //So we can work with i 
        if(added)
        {
           
            ans +=  1 ; 
        }
       }
        
          start_3_cli[index] =ans  ; 

      }
  }

      
//I have counted the number of 3 cliques , lets construct them 



__global__ void make_three_cliques (ll* edges ,ll*  two_cliques,ll    n ,ll  m , ll* start_3_cli , ll* three_cliques   ) 
      {
          
          ll index = blockDim.x * blockIdx.x + threadIdx.x ;
       if(index<m/2)
       {

           
           ll start_index = start_3_cli[index]  ; 

          ll ans = 0 ;    
          ll starting_element =  start_ele(two_cliques , m/2 , two_cliques[index*2+1]  ) ; 

          if(starting_element == -1 )
          return  ; 
          for(    ; starting_element < m/2 && two_cliques[starting_element*2] ==  two_cliques[index*2+1] ; starting_element++  )
          {
           ll i = two_cliques[starting_element*2+1] ;
    
           bool added = 1 ; 
          for(ll j = 0 ; j < 2; j++ )
          {
              if(check_edge(edges , 0 , m-1,  i , two_cliques[index*2+j] ) == 0 )
              {
                  added = 0 ;
                  break ;  
              }
          }
       //So we can work with i 
        if(added)
        {
                    three_cliques[start_index*3] =  two_cliques[index*2] ; 
                    three_cliques[start_index*3+1] =  two_cliques[index*2+1] ; 
                    three_cliques[start_index*3+2] =  i ; 
                    start_index++ ;  
            ans +=  1 ; 
        }
       }
        

        
        
        }

      }
//-----------------------------------------------------------------
//I will have to find the solution



__device__ ll find_starting_three(ll* three_cliques , ll left , ll right  ,  ll starting_element )
{
    ll ans = -1 ; 
 while(left<=right)
 {
     ll center = (right - left)/2 + left ; 
  if(three_cliques[center*3] == starting_element )
  {
      ans = center ; 
  }
  if(three_cliques [center*3] >= starting_element )
  {
      right = center -1 ; 
  } else 
  left = center +1 ; 
 }
 return ans ; 

}
__global__ void solution_greater (ll* edges ,ll*  two_cliques, ll* three_cliques,ll  n ,ll m ,ll k  ,ll three_size ,ll* update)
{
       ll index = blockDim.x * blockIdx.x + threadIdx.x ;

     if(index < three_size)
     {  
    //     printf("%d\n", index) ; 
       update[index]= 0 ; 
          if(k==4)
          {
              
            
            ll ans = 0 ;    
            ll starting_element =  start_ele(two_cliques , m/2 , three_cliques[index*3+2]  ) ; 

            if(starting_element == -1 )
            return  ; 
           
            for(    ; starting_element < m/2 && two_cliques[starting_element*2] ==  three_cliques[index*3+2] ; starting_element++  )
            {
             
             
              int counted = 1 ; 
                for(int j =0 ; j < 2 ; j++)
                {
                    for(int k = 0; k <3 ; k++)
                    {
                        //three clique .. two_clique 
                     if(j!=0 && k!= 2 )
                        if(check_edge(edges,0,m-1,two_cliques[starting_element*2+j ] , three_cliques[index*3+k])== 0   )
                        {
                            counted =0 ; 
                            break ;    
                        }
                    }
                }
             update[index] += counted ; 




            }
        



          }else if (k==5)
         {
             
             ll starting_element = find_starting_three(three_cliques , 0 , three_size -1 ,  three_cliques[index*3+2]) ; 
             if(starting_element == -1 )
             {
                 return ; 
             }
          
          for( ; starting_element< three_size && three_cliques[starting_element*3] ==  three_cliques[index*3+2] ; starting_element++ )
          {
                int counted = 1 ; 
                for(int j =0 ; j < 3 ; j++)
                {
                    for(int k = 0; k <3 ; k++)
                    {
                        //three clique .. two_clique 
                     if(j!=0 && k!=2  )
                        if(check_edge(edges,0,m-1, three_cliques[starting_element*3+j ] , three_cliques[index*3+k])== 0    )
                        {
                            counted =0 ; 
                            break ;    
                        }
                    }
                }
             update[index] += counted ; 

          }



         }else if(k==6)
       {
           
            
             
             ll starting_element = find_starting_three(three_cliques , 0 , three_size -1 ,  three_cliques[index*3+2]) ; 
             if(starting_element == -1 )
             {
                 return ; 
             }
          
          for( ; starting_element< three_size && three_cliques[starting_element*3] ==  three_cliques[index*3+2] ; starting_element++ )
          {
                int counted = 1 ; 
                for(int j =0 ; j < 3 ; j++)
                {
                    for(int k = 0; k <3 ; k++)
                    {
                        //three clique .. two_clique 
                     if(j!=0 && k!=2  )
                        if(check_edge(edges,0,m-1, three_cliques[starting_element*3+j ] , three_cliques[index*3+k])== 0    )
                        {
                            counted =0 ; 
                            break ;    
                        }
                    }
                }
           if(counted )
           {


                   LocalVector<ll> dum  ; 
                    dum.add( three_cliques[index*3]  ) ;
                    dum.add( three_cliques[index*3+1]  ) ;
                    dum.add( three_cliques[index*3+2]  ) ;
                    dum.add( three_cliques[starting_element*3+1]  ) ;
                    dum.add( three_cliques[starting_element*3+2]  ) ;
                    dum.add(0) ; 

              update[index] +=  sol(edges , dum , k , 5 , n , m , two_cliques ) ; 
         
 
  
           }
             
          }


       }else if(k==7)
       {
           

  
             ll starting_element = find_starting_three(three_cliques , 0 , three_size -1 ,  three_cliques[index*3+2]) ; 
             if(starting_element == -1 )
             {
                 return ; 
             }
          
          for( ; starting_element< three_size && three_cliques[starting_element*3] ==  three_cliques[index*3+2] ; starting_element++ )
          {
                int counted = 1 ; 
                for(int j =0 ; j < 3 ; j++)
                {
                    for(int k = 0; k <3 ; k++)
                    {
                        //three clique .. two_clique 
                     if(j!=0 && k!=2  )
                        if(check_edge(edges,0,m-1, three_cliques[starting_element*3+j ] , three_cliques[index*3+k])== 0    )
                        {
                            counted =0 ; 
                            break ;    
                        }
                    }
                }
           if(counted == 1 )
           {
               
                   LocalVector<ll> dum  ; 
                    dum.add( three_cliques[index*3]  ) ;
                    dum.add( three_cliques[index*3+1]  ) ;
                    dum.add( three_cliques[index*3+2]  ) ;
                    dum.add( three_cliques[starting_element*3+1]  ) ;
                    dum.add( three_cliques[starting_element*3+2]  ) ;
  
            
        
                  ll second_starting = find_starting_three(three_cliques , 0 , three_size -1 , dum[4] )  ;
                  if(second_starting == -1 )
                    continue ; 
                  for( ; second_starting < three_size && three_cliques[second_starting*3] == dum[4] ; second_starting++ )
                  {
                    bool okd = 1 ; 
                   for(ll gi = 0 ; gi <5 ;gi++)
                   {
                       for(ll gj  = 1  ; gj < 3 ; gj++ )

                       {
                           if(check_edge(edges , 0 , m-1, dum[gi] ,  three_cliques[second_starting*3+gj] )==0 )
                           okd = 0 ; 
                       }
                   }
                    update[index] += okd  ; 
                       
                  }
          }
          
          }



           
       }
      
     }
      
}


//-------------------------------------------------------------------------------------------------------------
//Values for 6 && 7 



__global__ void make_four_clique (ll* edges ,ll*  two_cliques, ll* three_cliques,ll  n ,ll m ,ll k  ,ll three_size ,ll* four_clique_sizes , ll* four_cliques)
{

     ll index = blockDim.x * blockIdx.x + threadIdx.x ;

     if(index < three_size)
     {  
          
              
            ll start_writing = four_clique_sizes[index] ; 
              
            ll ans = 0 ;    
            ll starting_element =  start_ele(two_cliques , m/2 , three_cliques[index*3+2]  ) ; 

            if(starting_element == -1 )
            return  ; 
           
            for(    ; starting_element < m/2 && two_cliques[starting_element*2] ==  three_cliques[index*3+2] ; starting_element++  )
            {
             
             
              int counted = 1 ; 
                for(int j =0 ; j < 2 ; j++)
                {
                    for(int k = 0; k <3 ; k++)
                    {
                        //three clique .. two_clique 
                     if(j!=0 && k!= 2 )
                        if(check_edge(edges,0,m-1,two_cliques[starting_element*2+j ] , three_cliques[index*3+k])== 0   )
                        {
                            counted =0 ; 
                            break ;    
                        }
                    }
                }
             if(counted)
             {
                 four_cliques[start_writing*4 ] = three_cliques[index*3] ; 
                 four_cliques[start_writing*4 +1] = three_cliques[index*3+1] ; 
                 four_cliques[start_writing*4 +2 ] = three_cliques[index*3+2] ; 
                 four_cliques[start_writing*4 +3 ] = two_cliques[starting_element*2+1] ; 
                start_writing++ ; 
             }
             

            }

      }


}

//------------------------------------------------------------------------

//Make  for 6 & 7 ------------------------------------------------------


__device__ ll find_starting_four(ll* four_cliques , ll left  ,ll right  ,  ll search )
{
    ll ans = -1 ; 
     while(left <= right )
     {
        ll center = (right - left)/2 + left ; 
        if(four_cliques[center*4] ==  search )
        ans = center ; 
        if(four_cliques[center*4] >= search)
        {
            right = center - 1; 
        }else 
      left = center +1 ;
     }
 return ans ; 

}
__global__  void solution_sizx_seven (ll* edges ,ll* four_cliques , ll* three_cliques ,ll  m ,ll number_of_four_cli  ,ll three_size ,ll*  update , ll k  ) 
{
    
     ll index = blockDim.x * blockIdx.x + threadIdx.x ;
    if(index < number_of_four_cli)
    {
        update[index] = 0 ; 
        if(k==6)
        {
            //4 and 3 







               ll starting_element = find_starting_three(three_cliques , 0 , three_size -1 ,  four_cliques[index*4+3])  ;  // last element  
                if(starting_element == -1 )
                {
                 return ; 
                }
          
          for( ; starting_element< three_size && three_cliques[starting_element*3] ==  four_cliques[index*4+3] ; starting_element++ )
          {
                int counted = 1 ; 
                for(int j =0 ; j < 3 ; j++)
                {
                    for(int k = 0; k <4 ; k++)
                    {
                        //four clique .. two_clique 
                     if(j!=0 && k!=3  )
                        if(check_edge(edges,0,m-1, three_cliques[starting_element*3+j ] ,four_cliques[index*4+k])== 0    )
                        {
                            counted =0 ; 
                            break ;    
                        }
                    }
                }
             update[index] += counted ; 

          }
        }else 
     {
               ll starting_element = find_starting_four(four_cliques , 0 , number_of_four_cli -1 ,  four_cliques[index*4+3])  ;  // last element  
                if(starting_element == -1 )
                {
                 return ; 
                }
          
          for( ; starting_element< number_of_four_cli && four_cliques[starting_element*4] ==  four_cliques[index*4+3] ; starting_element++ )
          {
                int counted = 1 ; 
                for(int j =0 ; j < 4 ; j++)
                {
                    for(int k = 0; k <4 ; k++)
                    {
                        //four .. three 
                     if(j!=0 && k!=3  )
                        if(check_edge(edges,0,m-1, four_cliques[starting_element*4+j ] ,four_cliques[index*4+k])== 0    )
                        {
                            counted =0 ; 
                            break ;    
                        }
                    }
                }
             update[index] += counted ; 

          }
         
     }
     
    }
}



//------------------------------------------------------------------

void  take_input( vector<ll> & edges , vector<ll> & two_clique , ll m )
{
  vector<pair<ll,ll>> unpacked_two_cli ; 
  for(ll i=0;i<m;i++)
  {
      
      ll a , b ; indata>>a>>b ;
      if(a>b)
      swap(a,b) ; 
      unpacked_two_cli.push_back({a,b}) ; 
  }
 
  sort(unpacked_two_cli.begin() , unpacked_two_cli.end()) ;
  vector<pair<ll,ll>> unpacked_edges = unpacked_two_cli; 
  for(auto i :  unpacked_two_cli )
  {
      unpacked_edges.push_back({i.second , i.first}) ; 
  }
  sort(unpacked_edges.begin(),  unpacked_edges.end()); 
 for(ll i=0;i<m;i++)
 {
     
 two_clique[i*2] = unpacked_two_cli[i].first ;
  
 two_clique[i*2+1] = unpacked_two_cli[i].second ;
 }


 for(ll i=0;i<2*m;i++)
 {
     
 edges[i*2] = unpacked_edges[i].first ;
  
 edges[i*2+1] = unpacked_edges[i].second ;
 }




}



//--------------------------------------------------------------------------------
ll solution(string file_name,  ll k )
{   
    
    indata.open(file_name); // opens the file
    if(!indata) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
    }
    
 //   freopen("com-youtube.ungraph.txt","r",stdin); //as20000102
    
    ll n = 1 ;
    ll  m ; indata>> m ; 

    vector<ll> edges(4*m) , two_cliques(2*m) ;
    take_input(edges ,  two_cliques , m ) ; 
    for(auto  i : two_cliques)
    n = max(i , n ) ; 
      if(k==2)
      return m ; 
 

    n++  ;
    m*=2 ; 
 
if(k==1)
return n-1 ; 
 if(k==2)
 return m/2 ; 
 else 
    {
      
      /*--------------------------------------------------------------Copy the edges and two cliques----------------------------------------------*/
      ll *device_edges  , *device_updates , *device_two_cliques ;
  

      cudaMalloc( &device_edges ,  sizeof(ll)*m*2)  ;
      cudaMalloc( &device_two_cliques , sizeof(ll)*(m) )  ;

      cudaMemcpy( device_edges , edges.data() , sizeof(ll)*m*2 , cudaMemcpyHostToDevice )  ;   
      cudaMemcpy(  device_two_cliques , two_cliques.data() , sizeof(ll)*m , cudaMemcpyHostToDevice )  ;   

      /*------------------------------------------------------------------------------------------------------------------------------------------*/
      






      /*---------------------------------------------------------------Count number of three cliques----------------------------------------------*/
      vector<ll> start_3_cli (m/2)  ;
      ll * device_start_3_cli ; 
      cudaMalloc( &device_start_3_cli , sizeof(ll)*m/2) ; 
      
      
      // Invoke kernel
      int threadsPerBlock = 100 ;
      ll blocksPerGrid = ( (m/2 )+ threadsPerBlock - 1) / threadsPerBlock;
      count_three_cliques<<<blocksPerGrid, threadsPerBlock>>> (device_edges , device_two_cliques,   n , m , device_start_3_cli  );
       

      cudaMemcpy(start_3_cli.data() , device_start_3_cli, sizeof(ll)*(m/2) , cudaMemcpyDeviceToHost);
      ll current_size =0 ; 
     
      for(ll i= 0;i<m/2 ;i++)
      {    
         ll temp =start_3_cli[i]  ; 
         start_3_cli[i] = current_size ; 
        current_size += temp ; 
      }
      /*----------------------------------------------------------------------------------------------------------------------------------------*/
     
     
      if(k==3)
      {
          cudaFree(device_start_3_cli) ; 
          cudaFree(device_edges) ; 
          cudaFree(device_updates) ;
          cudaFree(device_two_cliques); 
          return current_size ; 
      }


     /*-------------------------------------------------Copy all the three cliques-------------------------------------------------------------*/

      if(current_size ==0 )
        return 0 ; 
      ll* three_cliques  ; 
      cudaMalloc( &three_cliques , sizeof(ll)*current_size*3) ;
       
      cudaMemcpy(device_start_3_cli ,start_3_cli.data()  , sizeof(ll)*(m/2) , cudaMemcpyHostToDevice )  ;   
    
       make_three_cliques<<<blocksPerGrid, threadsPerBlock>>> (device_edges , device_two_cliques,   n , m , device_start_3_cli , three_cliques   );

    /*-----------------------------------------------------------------------------------------------------------------------------------------*/
     
    ll sol_ans =0 ;
 

    if(k <6)
    {
     /*----------------------------------------------------------------Find the solutions-------------------------------------------------------*/
        threadsPerBlock = 1024 ;
       blocksPerGrid = ( (current_size )+ threadsPerBlock - 1) / threadsPerBlock;
      
      
           vector<ll>updates(current_size) ;
         cudaMalloc( &device_updates , sizeof(ll)*(current_size) )  ;

     
        solution_greater<<<blocksPerGrid, threadsPerBlock>>> (device_edges , device_two_cliques, three_cliques,    n , m , k ,current_size , device_updates );
         cudaMemcpy(updates.data() , device_updates, sizeof(ll)*(current_size) , cudaMemcpyDeviceToHost);

      for(auto i :  updates) 
        sol_ans +=  i ;
     
    /*------------------------------------------------------------------------------------------------------------------------------------------*/
        
    }else 
     {
         
         /*----------------------------------------------------------------Find the solutions-------------------------------------------------------*/
        threadsPerBlock = 1024 ;
       blocksPerGrid = ( (current_size )+ threadsPerBlock - 1) / threadsPerBlock;
      
        ll *four_clique_size ; 
        cudaMalloc( &four_clique_size , sizeof(ll)*(current_size) )  ;

     
        solution_greater<<<blocksPerGrid, threadsPerBlock>>> (device_edges , device_two_cliques, three_cliques,    n , m , 4 ,current_size , four_clique_size ); // We only need 4 cliques 
        vector<ll> four_cli_si(current_size ) ; 
        cudaMemcpy(four_cli_si.data() , four_clique_size, sizeof(ll)*(current_size) , cudaMemcpyDeviceToHost);
  
        ll number_of_four_cli =  0 ;
        for(ll i= 0;i<current_size;i++)
        {
            ll temp = four_cli_si[i] ; 
            four_cli_si[i] = number_of_four_cli ; 
            number_of_four_cli += temp ;
        }
      
        cudaMemcpy(  four_clique_size,four_cli_si.data() ,  sizeof(ll)*(current_size) , cudaMemcpyHostToDevice) ;
        
        //sol_ans = number_of_four_cli ; 
            /*-------------------------------------------------------------Now make the four clique----------------------------------------------*/
              ll* four_cliques ; 
              cudaMalloc(&four_cliques, sizeof(ll)*4*number_of_four_cli ) ; 
              make_four_clique<<<blocksPerGrid, threadsPerBlock>>> (device_edges , device_two_cliques, three_cliques,    n , m , k ,current_size , four_clique_size , four_cliques );
            /*------------------------------------------------------------------------------------------------------------------------------------*/
      


            /*----------------------------------------------------------Now ask the values of 6 and 7 cliques-------------------------------------*/
            ll *device_updates ; 
              blocksPerGrid = ( (number_of_four_cli )+ threadsPerBlock - 1) / threadsPerBlock;
              cudaMalloc(&device_updates , sizeof(ll)*number_of_four_cli ) ;
          
                        vector<ll> arsd(number_of_four_cli) ; 
                       solution_sizx_seven<<<blocksPerGrid,threadsPerBlock>>> (device_edges , four_cliques , three_cliques , m , number_of_four_cli  , current_size,  device_updates , k ) ; 
                    
                      cudaMemcpy(arsd.data() , device_updates, sizeof(ll)*number_of_four_cli , cudaMemcpyDeviceToHost);  

                      
                          for(ll i =0;i<  number_of_four_cli ;i++)
                              sol_ans += arsd[i] ;                   
      
                      

            /*------------------------------------------------------------------------------------------------------------------------------------*/



      
    /*------------------------------------------------------------------------------------------------------------------------------------------*/
     


     }
     

      cudaFree(device_updates);
      cudaFree(device_edges);
      cudaFree(device_two_cliques);
      cudaFree(three_cliques);
      cudaFree(device_start_3_cli);
     
      return sol_ans ; 
     
    }
}



int main()
{

      string file_name  ; cin>>file_name ; 
    
    ll k ; cin>>k  ;
    time_t x1, x2 ;
    x1 = time(NULL);
    
    string sd =""  ;
    ll ans = solution(file_name, k) ;
    x2 = time(NULL) ; 
    cout<<"\n"<<x2-x1<<" " << k<<" "<<ans  ; 

 

}