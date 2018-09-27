// Copyright (C) 2004, 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: MyNLP.hpp 511 2005-08-26 22:20:20Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-11-05

#ifndef __MYNLP_HPP__
#define __MYNLP_HPP__

#define DEBUG_hd 0
//#include "IpTNLP.hpp"

//using namespace Ipopt;

#include "placedb.h"
#include "threadpool.h"
#include <set>
#include <vector>
#include <string>
#include <string.h>
using namespace std;


class MyNLP 
{

public:
  MyNLP(){}
  MyNLP( CPlaceDB& );

  virtual ~MyNLP();

  bool get_nlp_info(int& n, int& m, int& nnz_jac_g, int& nnz_h_lag );
  bool get_bounds_info(int n, double* x_l, double* x_u, int m, double* g_l, double* g_u);
  bool get_starting_point(int n, bool init_x, double* x, bool init_z, double* z_L, double* z_U,
                          int m, bool init_lambda, double* lambda);
  bool eval_f(int n, const double* x, const double* expX, bool new_x, double& obj_value);
  bool eval_f_HPWL(int n, const double* x, const double* expX, bool new_x, double& obj_value);
  bool eval_grad_f( int n, const double* x, const double* expX, bool new_x, double* grad_f);
  
  // return true is placement is legal
  bool MySolve( double, double target_density, int currentLevel, bool noRelaxSmooth,threadpool_t **pool = nullptr );	// solver setup
  
  int _potentialGridR;
  int m_potentialGridSize;
  double m_targetUtil;
  double target_nnb;

  bool   m_lookAheadLegalization;
  bool   m_earlyStop;
  bool   m_topLevel;
  int    m_smoothR;
  bool   m_lastNLP;
  bool   m_useBellPotentialForPreplaced;
  double m_smoothDelta;

  double* x;
  
private:
  bool GoSolve( double, double target_density, int currentLevel );	// real solver
  
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually 
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *  
   */
  //@{
  //  MyNLP();
  MyNLP(const MyNLP&); 
  MyNLP& operator=(const MyNLP&);
  //@}


  CPlaceDB* m_pDB;
  vector< pair<int,int> > _cellPair;
  void UpdateBlockPosition( const double* x );
  void BoundX( const int& n, double* x, double* x_l, double* x_h );
  inline void BoundX( const int& n, double* x, double* x_l, double* x_h, const int& i );
  void LineSearch( const int& n, /*const*/ double* x, double* grad_f, double& stepSize );
  void AdjustForce( const int& n, const double* x, double* grad_f );
  void AdjustForce( const int& n, const double* x, vector<double> grad_wl, vector<double> grad_potential );
  void FindBeta( const int& n, const double* grad_f, const double* last_grad_f, double& beta );
  double _weightWire;
  double _weightDensity;

  double* xBest;   // look ahead legalization 
  double* x_l;	   // lower bound
  double* x_u;     // upper bound
  double* _expX;   // exp(x)
  double* _expPins;
public:
  double* xBak;
  double* xBak2;
private:
  vector<double> grad_wire;
  vector<double> grad_potential;
  int m_ite;
  double m_currentStep;
  
public:
  double m_incFactor;
  double m_weightWire;

private:
  // wirelength related functions
  void calc_sum_exp_using_pin(
          const vector<int>::const_iterator& begin, const vector<int>::const_iterator& end,
          const double* x, const double* expX,
          double& sum_exp_xi_over_alpha, double& sum_exp_inv_xi_over_alpha,
          double& sum_exp_yi_over_alpha, double& sum_exp_inv_yi_over_alpha, int id=-1 );
  double GetWL( const int& n, const double* x, const double* expX, const double& alpha );
  void   UpdateExpValueForEachCell( const int& n, const double* x, double* expX, const double& alpha );
  void   UpdateExpValueForEachPin( const int& n, const double* x, double* expPins, const double& alpha );
  void   UpdateNetsSumExp( const double* x, const double* expX );

  double m_posScale;
  double _alpha;
  vector<bool> m_usePin;
  void SetUsePin(); 
  void InitModuleNetPinId();
  vector< vector<int> > m_moduleNetPinId;

  vector<double> m_nets_sum_exp_xi_over_alpha;
  vector<double> m_nets_sum_exp_yi_over_alpha;
  vector<double> m_nets_sum_exp_inv_xi_over_alpha;
  vector<double> m_nets_sum_exp_inv_yi_over_alpha;

  vector<double> m_nets_sum_p_x_pos;
  vector<double> m_nets_sum_p_y_pos;
  vector<double> m_nets_sum_p_inv_x_pos;
  vector<double> m_nets_sum_p_inv_y_pos;
  vector<double> m_nets_sum_p_x_neg;
  vector<double> m_nets_sum_p_y_neg;
  vector<double> m_nets_sum_p_inv_x_neg;
  vector<double> m_nets_sum_p_inv_y_neg;
  
  // diffusion bin
  vector< vector<double> > m_binForceX;
  vector< vector<double> > m_binForceY;
  void UpdateBinForce();
  void GetDiffusionGrad( const double* x, const int& i, double& gradX, double& gradY );
  
  // potential grid related variables/functions
  double m_alwaysOverPotential;
  vector< double > _cellPotentialNorm;		// cell potential normalization factor
  vector< vector<double> > m_gridPotential;
  double m_potentialGridWidth;
  double m_potentialGridHeight;
  double _potentialRY;
  double _potentialRX;
  //double _expPotential;
  double GetNonZeroGridPercent();
  double GetMaxPotential();
  double GetAvgPotential();
  double GetTotalOverPotential();   // 2006-03-01
  void   OutputPotentialGrid( string filename );	// for gnuplot
  void   OutputDensityGrid( string filename );		// for gnuplot
  void   PrintPotentialGrid();
  void   GetPotentialGrad( const double* x, const int& i, double& gradX, double& gradY );
  void   CreatePotentialGrid();
  void   ClearPotentialGrid();
  void   UpdatePotentialGrid( const double* x );
  void   UpdatePotentialGridBase( const double* x );	    // compute preplaced block potential
  void   SmoothPotentialBase( const double& delta );	    // 2006-03-04
  double GetGridWidth();
  void   GetGridCenter( const int& gx, const int& gy, double& x1, double& y1 );
  double GetXGrid( const int& gx );
  double GetYGrid( const int& gy );
  double GetDensityPanelty();
  double GetPotentialToGrid( const double& x1, const int& gx );         // 1D
  double GetGradPotentialToGrid( const double& x1, const int& gx );	// 1D
  void   GetClosestGrid( const double& x1, const double& y1, int& gx, int& gy );
  struct potentialStruct    
  {
      potentialStruct( const int& x, const int& y, const double& p ) 
	  : gx(x), gy(y), potential(p)
	  {}
      int gx;
      int gy;
      double potential;
  };
  void UpdateExpBinPotential( double utl );	// 2006-03-14
  vector< vector<double> > m_basePotential;	// 2006-03-03 (donnie) preplaced block potential 
  vector< vector<double> > m_binFreeSpace;	// 2006-03-16 (donnie) free space in the bin 
  vector< vector<double> > m_basePotentialOri;	// 2006-03-03 (donnie) preplaced block potential 
  vector< vector<double> > m_expBinPotential;	// 2006-03-14 (donnie) preplaced block potential 
  
  
  // bell-shaped functions
  inline double GetPotential( const double& x1, const double& x2, const double& r );
  inline double GetPotential( const double& x1, const double& x2, const double& r, const double& w );
  inline double GetGradPotential( const double& x1, const double& x2, const double& r );
  inline double GetGradPotential( const double& x1, const double& x2, const double& r, const double& w );
  inline double GetGradGradPotential( const double& x1, const double& x2, const double& r ); 

  
  // density grid related functions
  vector< vector<double> > m_gridDensity;
  vector< vector<double> > m_gridDensitySpace;	// avaliable space in the bin
  double m_alwaysOverArea;
  //double m_totalMovableModuleArea;
  //double m_totalFixedModuleArea;
  double m_totalFreeSpace;
  double m_gridDensityWidth;
  double m_gridDensityHeight;
  double m_gridDensityTarget;
  void   UpdateDensityGrid( const int& n, const double* x );
  void   UpdateDensityGridSpace( const int& n, const double* x );
  void   CheckDensityGrid();
  void   CreateDensityGrid( int nGrid );
  void   ClearDensityGrid();
  double GetDensityGridPanelty();
  double GetMaxDensity();
  double GetAvgOverDensity();
  double GetTotalOverDensity();
  double GetTotalOverDensityLB();
  double GetNonZeroDensityGridPercent();
  void   GetDensityGrad( const double* x, const int& i, double& gradX, double& gradY ); // UNUSE

 /*************************************************flpeng************************************************/
public:
  vector<double> sub_grad_wire;
  vector<double> sub_grad_potential;


  double sg_CalcHPWL(const double* x);
  void sg_findBeta(const int& n, const double* grad_f, const double* last_grad_f,const double* real_last_grad_f, double& beta );
  double getRandNum(double x1, double x2);
  bool sg_eval_f(const double* x, double& objValue);
  void sg_eval_grad_f(int n, const double* x, double* sub_grad_f);
  //simplfied method  to update the subgradient of wirelength, only calculate the max or min module
  void sg_updateSubGradWire_net2_s(const int &netID, const double *x,bool isPart = false);
  void sg_updateSubGradWire_netBig_s(const int &netID, const double *x,bool isPart = false);
  void sg_getFirstModuleSubWireGrad(const int& id1, const int& id2, const double& x1, const double& x2, const double& y1, const double& y2,bool isPart = false);
  void sg_getXModuleSubWireGrad(const int& id1, const int& id2, const double& x1, const double& x2,bool isPart );
  void sg_getYModuleSubWireGrad(const int& id1, const int& id2, const double& y1, const double& y2,bool isPart);

  /************************************************flpeng************************************************/


  /********************************************hdgao added**********************************************/
    //serialization conjucated gradient
public:
  vector< double > part_grad_wire;
  vector< double > part_grad_potential;
  vector< double > part_grad;

  vector< vector<int> > nets_cluster;
  vector< vector<int> > nets_black_modules;
  vector< vector<int> > nets_cluster_modules_bounds;

  vector< vector<int> > nets_cluster_WL;
  vector< vector< vector<int> > > nets_cluster_Den;
  vector<int> modules;// for parallel
  vector<bool> isread;


  vector< int > nets_random;
  vector< int > modules_random;



  static void CG(MyNLP* cur_mynlp,int size_x, double *x_bak, bool is_2D);
  static void CG(MyNLP* cur_mynlp,int size_x, double *x_bak, bool is_2D,// process some specified nets cluster of nets_cluster
                const vector<int>& specified_nets_cluster);
  static void CG(MyNLP* cur_mynlp,int size_x, double *x_bak, bool is_2D,// process some specified nets cluster of nets_cluster
        const vector< vector<int> >& cur_nets_cluster,
        const vector< vector<int> >& cur_nets_black_modules,
              vector< vector<int> > cur_nets_cluster_modules_bounds,
          const vector<int>& specified_nets_cluster);

  static void GD(MyNLP* cur_mynlp,int size_x, double *x_bak, bool is_2D,// process some specified nets cluster of nets_cluster
        const vector< vector<int> >& cur_nets_cluster,
        const vector< vector<int> >& cur_nets_black_modules,
              vector< vector<int> > cur_nets_cluster_modules_bounds,
          const vector<int>& specified_nets_cluster);

  static void task_CG(void* args);
  static void task_GD(void* args);
  friend void shell_task_CG(void(*fp)(void*),void* args);
  void MyNLP_copy( MyNLP* source);
  void MyNLPCopy(CPlaceDB&);
  double PL_step;
// parallel proceeding
public:
  static void SGD_Density(MyNLP* cur_mynlp,size_t size_x,const vector<int>& specified_nets_cluster,
                  const vector<int>& bin,double stepsize,double maxloop);
  static void SGD_WireLength(MyNLP* cur_mynlp,size_t size_x,const vector<int>& specified_nets_cluster,
                  double stepsize,double maxloop);
  static void task_Density(void* args);
  static void task_WireLength(void* args);

  //density
  void DensityGrad(vector<double>& grad,const vector<int> netsId,const vector<int>& bin = vector<int>());
  void DensityGradByModule(vector<double>& grad,const vector<int>& modules,const vector<int>& bin = vector<int>());
  void DensityAdjustForce(size_t size_x,vector<double>& grad,const vector<int> netsID);
  void DensityBoundX(size_t size_x);
  void DensityUpdatePotentialGrid(const vector<int> netsID,vector<int> bin);
  void DensityUpdatePotentialGridByModule(const vector<int>& modules,const vector<int>& bin = vector<int>());
//  void DensityUpdateExpValueForCell(const vector<int> netsID,vector<int> bin);
//  void DensityUpdateExpValueForPin(const vector<int> netsID);
//  void DensityUpdateNetsSumExp(const vector<int> netsID);
  // wirelength
  void WireLengthGrad(vector<double>& grad,const vector<int> netsID);
  void WireLengthGradByModule(vector<double>& grad,const vector<int>& modules);
  void WireLengthAdjustForce(size_t size_x,vector<double>& grad,const vector<int> netsID);
  void WireLengthBoundX(size_t size_x);
//  void WireLengthUpdatePotentialGrid();
  void WireLengthUpdateExpValueForCell(const vector<int> netsID,size_t size_x);
  void WireLengthUpdateExpValueForCellByModule(const vector<int> modules,size_t size_x);
  void WireLengthUpdateExpValueForPin(const vector<int> netsID);
  void WireLengthUpdateNetsSumExp(const vector<int> netsID);
  // helper
  void GetGrid(int bin_left,int bin_bottom,double bucket_width,double bucket_height,vector<int>& bin);
  void GetDensityGrad(const size_t moduleId,const vector<int> bin,double &gradx,double &grady);
  void DensityUpdateOneModulePotentialGrid(size_t moduleId,vector<int> bin);
  void GetModulesByNets(const vector<int> netsID,vector<int>& modules);
  void UpdateModulePotentialGridByDiff(size_t moduleId,vector<int> bin);



private:
  bool SLSolve_bak(double wWire, double target_density, int currentLevel);
  bool SLSolve_smooth_bak(double wWire, double target_density, int currentLevel);
  bool SLSolve(double wWire, double target_density, int currentLevel);
  bool SLSolve_nets_cluster(double wWire, double target_density, int currentLevel);
  // serialization process
  bool SLSolve_smooth(double wWire, double target_density, int currentLevel);
  // Parallel process
  bool PLSolve_smooth(double wWire, double target_density, int currentLevel,threadpool_t* pool);
  static void dummy_task(void *arg);
  void CG_mynlp(int size_x,double* x_bak,bool is_2D);
  void CG_mynlp(int size_x,double* x_bak,bool is_2D,const vector<int>& specified_nets_cluster);

  // hdgao added 15/09
  bool ParallelSolve(double wWire, double target_density, int currentLevel,threadpool_t** pool);



  void part_eval_grad_f(const double* x, double* sub_grad_f, int i_net, int num_nets = 1);
  double part_CalcHPWL(const double* x ,int i_nets);
  bool  test_part_eval_grad_f(const double* x);

  bool smooth_part_eval_grad_f(const double* x,double* part_grad_f,int i_nets,int num_nets);
  bool smooth_part_eval_grad_f(const double* x,
                               const vector<int> cur_nets_cluster_modules_bounds, const vector<int> cur_black_module,
                               vector<int>& RegionModule,
                               vector<double> &grad_f, int flag = true);

  void part_UpdatePotentialGrid(const double *x,const double *x_bak,const vector<int>& moduleId,const int num_module);
  void OneModulePotentialGrid(const int cur_moduleId,const double* position,vector< vector<double> >& potential_record);
  void part_UpdatePotentialGrid(const double* x);

  bool nets_cluster_analysics();
  void getBlackModule(int nets_clusterId);
  void part_eval_grad_f(const double* x, int net_cluster_id, vector<int>& RegionModule,
                        vector<double> &grad_f, int flag = true);
 // void add_region_module(const vector<int>& WMIds,vector<int>& RMIds);
  void updateBlackBounds(int net_cluserId);
  void updateBlackBounds(int net_cluserId,
                         const vector< vector<int> >& cur_nets_cluster,
                         const vector< vector<int> >& cur_nets_black_modules,
                         vector< vector<int> >& cur_nets_cluster_modules_bounds);

  void UpdateOneModulePotentialGrid(int moduleId,vector<int>& RegionModule,int region_gx1,int region_gx2,int region_gy1,int region_gy2,
                                    double& partGradDensityX,double& partGradDensityY);
  void UpdateOneModulePotentialGrid_adjust(int moduleId,
                                           //vector<int>& RegionModule,
                                           int region_gx1,int region_gx2,int region_gy1,int region_gy2,
                                           //int n_RegionModule,
                                           double& partGradDensityX,double& partGradDensityY);
  void UpdateOneModulePotentialGrid(int moduleId,double& partGradDensityX,double& partGradDensityY);
  void test(threadpool_t *pool = nullptr);

  void PartAdjustForce(vector<double> &grad_f, const vector<int> &RegionModule);
  void PartAdjustForce_2D(vector<double> &grad_f, const vector<int> &RegionModule);
  void PartFindBeta(const vector<double>& grad_f,const vector<double>& last_grad_f,
                    const vector<int>& RegionModule,double& beta);
  void PartFindBeta_2D(const vector<double>& grad_f,const vector<double>& last_grad_f,
                    const vector<int>& RegionModule,double& beta1,double& beta2);
  void PartLineSearch(const vector<double>& grad_f,const vector<int>& RegionModule,double& stepsize);
  void PartLineSearch_2D(const vector<double>& grad_f,const vector<int>& RegionModule,double& stepsize1,double& stepsize2);

  void AdjustForce( const int& n, const double* x, vector<double>& grad_f );
  void AdjustForce_2D( const int& n, const double* x, vector<double>& grad_f );
  void LineSearch( const int& n, /*const*/ double* x, vector<double>& grad_f, double& stepSize );
  bool NCA(int moduleId,vector<bool>& is_nets_record,vector<int>& module_net);
  bool NCA_deeper();
  void get_center_position(double *x, double& c_x, double& c_y, const size_t netId, int flag);
  std::vector<int> randVector(size_t num);
  /********************************************hdgao added**********************************************/
  
};
/*
class stuffs_4_parallel{
public:
    stuffs_4_parallel();
    stuffs_4_parallel(int nets_size, int modules_size, int pins_size);
    vector<double> m_nets_sum_exp_xi_over_alpha;
    vector<double> m_nets_sum_exp_yi_over_alpha;
    vector<double> m_nets_sum_exp_inv_xi_over_alpha;
    vector<double> m_nets_sum_exp_inv_yi_over_alpha;

    vector<double> m_nets_sum_p_x_pos;
    vector<double> m_nets_sum_p_y_pos;
    vector<double> m_nets_sum_p_inv_x_pos;
    vector<double> m_nets_sum_p_inv_y_pos;
    vector<double> m_nets_sum_p_x_neg;
    vector<double> m_nets_sum_p_y_neg;
    vector<double> m_nets_sum_p_inv_x_neg;
    vector<double> m_nets_sum_p_inv_y_neg;

    vector<double> _cellPotentialNorm;

    double* _expX;   // exp(x)
    double* _expPins;
    double* x;
    double* x_l;
    double* x_u;



};
*/
/*
void CG(MyNLP* cur_mynlp,int size_x, double *x_bak, bool is_2D);
void CG(MyNLP* cur_mynlp,int size_x, double *x_bak, bool is_2D,// process some specified nets cluster of nets_cluster
//        const vector< vector<int> >& cur_nets_cluster,
//        const vector< vector<int> >& cur_nets_black_modules,
//        const vector< vector<int> >& cur_nets_cluster_modules_bounds,
        const vector<int>& specified_nets_cluster);
void task_CG(void* args);
*/
void shell_task_CG(void(*fp)(void*),void* args);

class args_class{
public:
    args_class(){}
    args_class(const args_class&);           // not implement
    args_class& operator=(const args_class&);// not implement

    args_class(MyNLP* source,
               const vector<int>& snc_Id,
               const unsigned int _size_x,
               const vector< vector<int> > _nets_cluster,
               const vector< vector<int> > _nets_black_modules,
               const vector< vector<int> > _nets_cluster_modules_bounds,
               bool _is_2D = false){

        //args_mynlp = new MyNLP();
        args_mynlp = source;

        if(x_bak == nullptr)
            x_bak = new double[_size_x];
        for(int i = 0;i < _size_x;i++){
            x_bak[i] = 0;
            args_mynlp->x[i] = args_mynlp->x[i];
            x_bak[i] = args_mynlp->x[i];
        }
        if(_nets_cluster.size() != 0){

            args_mynlp->nets_cluster.clear();
            args_mynlp->nets_black_modules.clear();
            args_mynlp->nets_cluster_modules_bounds.clear();

            args_mynlp->nets_cluster.resize(_nets_cluster.size());
            args_mynlp->nets_black_modules.resize(_nets_black_modules.size());
            args_mynlp->nets_cluster_modules_bounds.resize(_nets_cluster_modules_bounds.size());

            for(size_t i = 0;i < _nets_cluster.size();i++)
                args_mynlp->nets_cluster[i].assign(_nets_cluster[i].begin(),
                                       _nets_cluster[i].end());


            for(size_t i = 0;i < _nets_black_modules.size();i++)
                args_mynlp->nets_black_modules[i].assign(_nets_black_modules[i].begin(),
                                             _nets_black_modules[i].end());

            for(size_t i = 0;i < _nets_cluster_modules_bounds.size();i++)
                args_mynlp->nets_cluster_modules_bounds[i].assign( _nets_cluster_modules_bounds[i].begin(),
                                                       _nets_cluster_modules_bounds[i].end() );

        }


        size_x = _size_x;
        is_2D = _is_2D;

        args_snc_Id.assign(snc_Id.begin(),snc_Id.end());

    }
    args_class(MyNLP* source,const vector<int>& snc_Id,const unsigned int _size_x){

        //args_mynlp = new MyNLP();
        args_mynlp = source;
        size_x = _size_x;
        args_snc_Id.assign(snc_Id.begin(),snc_Id.end());
    }
        MyNLP* args_mynlp;
    vector<int> args_snc_Id;

    size_t size_x;
    unsigned int thread_number;
    double* x_bak;
    bool is_2D;
    double temp;
    double stepsize;
    size_t maxloop;
    vector<int> bin;
};


#endif
