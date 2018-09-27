#include <cmath>
#include <set>
#include <vector>
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <fstream>
#include <cstring>
#include <limits.h>
#include <unistd.h>
#include <stdlib.h>


using namespace std;

#include "placedb.h"
#include "MyNLP.h"
#include "smooth.h"
#include "TetrisLegal.h"
#include "placebin.h"
#include "ParamPlacement.h"

//double stepScaling = 0.2;
double truncationFactor = 1.0;
//double truncationFactor = 1.25;
//double truncationFactor = 0.5;

double time_grad_f = 0;
double time_f = 0;
double density;
double totalWL;

double time_wl = 0;
double time_update_grid = 0;
double time_grad_wl = 0;
double time_grad_potential = 0;
double time_update_bin = 0;
double time_up_potential = 0;
double time_exp = 0;
double time_sum_exp = 0;
double time_log = 0;

// for Euler method
double gradWL;
double stepSize;   


//for the precision  added by flpeng
const double sgEps = 1e-60;

//for parallel

pthread_mutex_t lock;
pthread_cond_t  cond_lock;
unsigned int job_done;
bool done_flag;
bool parallel_flag = true;
vector< vector< double > > save_x_threads;
double *Xmerge;



/*
void MyNLP::CreateExpTable()
{
    _expTablePrecision = 20000000;            // 20M * 8byte = 160M table
    _expTable.resize( _expTablePrecision+1 );

    double minValue = 1e10;
    double maxValue = -1e10;
    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
    if( m_pDB->m_modules[i].m_cx > maxValue )     maxValue =  m_pDB->m_modules[i].m_cx;
    if( m_pDB->m_modules[i].m_cy > maxValue )     maxValue =  m_pDB->m_modules[i].m_cy;
    if( m_pDB->m_modules[i].m_cx < minValue )     minValue =  m_pDB->m_modules[i].m_cx;
    if( m_pDB->m_modules[i].m_cy < minValue )     minValue =  m_pDB->m_modules[i].m_cy;
    }
    if( m_pDB->m_coreRgn.left   < minValue )    minValue = m_pDB->m_coreRgn.left;
    if( m_pDB->m_coreRgn.bottom < minValue )    minValue = m_pDB->m_coreRgn.bottom;
    if( m_pDB->m_coreRgn.top    > maxValue )    maxValue = m_pDB->m_coreRgn.top;
    if( m_pDB->m_coreRgn.right  > maxValue )    maxValue = m_pDB->m_coreRgn.right;

    maxValue /= _alpha;
    minValue /= _alpha;

    
    // TODO: check boundary
//    if( minValue < 0.0 )
//        minValue = 0.0;

    // TODO: maxValue < 0
    
    double step = ( maxValue - minValue ) / (_expTablePrecision);
    printf( " exp table: min %f(*alpha=%f)  max %f(*alpha=%f)  step %f\n",
        minValue, minValue*_alpha,
        maxValue, maxValue*_alpha,
        step );
    for( int i=0; i<_expTablePrecision+1; i++ )
    {
    _expTable[i] = exp( minValue + step * i );
    }
    _expTableMin = minValue;
    _expTableMax = maxValue;
    _expTableStep = step;
    //printf( "table[100] = %.10f\n", _expTable[100] );
    //printf( "exp(minValue+100*step) = %.10f\n", exp(minValue+step*100 ) );

}
*/

/*
double MyNLP::expTable( const double& value )
{
    double t_start = seconds();

    // TEST PRECISION && RUNTIME
    double res = exp(value);
    time_exp += seconds() - t_start;
    return res;



    
//    double absValue = fabs( value );
    double absValue = value;
    
    // BOUNDARY CHECK START =====================
    if( absValue < _expTableMin )
    {
    printf( "ERROR! value= %f (%f) min= %f\n", value, value*_alpha, _expTableMin );
    }
    if( absValue > _expTableMax )
    {
//	printf( "WARNING! value= %f (%f) max= %f\n", value, value*_alpha, this->_expTableMax );
    return exp( value );
    }
    //assert( absValue >= _expTableMin );
    //assert( absValue <= _expTableMax );
    // BOUNDARY CHECK END =======================
    
    
    // BACK LOOKUP
    //int index1 = static_cast<int>( ceil( (value - _expTableMin) / _expTableStep ) );
    //return _expTable[index1];		// 1% err

    // INTERPOLATION
    int index1 = static_cast<int>( floor( (absValue - _expTableMin) / _expTableStep ) );
    int index2 = index1 + 1 ;
    //assert( value >= _expTableMin + index1 * _expTableStep );
    //assert( value <  _expTableMin + index2 * _expTableStep );

    double diff = absValue - ( _expTableMin + index1 * _expTableStep );
    //assert( diff >= 0 );
    //assert( diff - _expTableStep <= 0 );

    double pos = _expTable[index1] +
       (_expTable[index2] - _expTable[index1]) * (diff) / _expTableStep;

    time_exp += seconds() - t_start;
//    if( value > 0 )
    return pos;
//    else
//	return 1.0 / pos;
}
*/

/*
void MyNLP::CheckExpTablePrecision()
{
    double maxValue = 1500 / _alpha;
    double step = 0.0000001;
    double diff = 0.0;
    for( double i=-maxValue; i<maxValue; i+=step )
    {
    double e1 = expTable( i );
    double e2 = exp(i);
    //printf( "e2 = %f   e1 = %f   diff = %f\n", e2, e1, e2-e1 );
    diff += (e2-e1)*(e2-e1);
    }
    printf( "CheckExpTablePrecision %f-%f %f err= %g\n", -maxValue, maxValue, step, diff );
    
    double t1 = seconds();
    double e1;
    for( double i=-maxValue; i<maxValue; i+=step )
    e1 = expTable( i );
    t1 = seconds() - t1;
    double t2 = seconds();
    for( double i=-maxValue; i<maxValue; i+=step )
    e1 = exp( i );
    t2 = seconds() - t2;
    printf( "time_compare = %f %f\n", t1, t2 );
}
*/
/*
stuffs_4_parallel::stuffs_4_parallel(int nets_size,int modules_size,int pins_size){
    m_nets_sum_exp_xi_over_alpha.resize( nets_size, 0 );
    m_nets_sum_exp_yi_over_alpha.resize( nets_size, 0 );
    m_nets_sum_exp_inv_xi_over_alpha.resize( nets_size, 0 );
    m_nets_sum_exp_inv_yi_over_alpha.resize( nets_size, 0 );

    m_nets_sum_p_x_pos.resize( nets_size, 0 );
    m_nets_sum_p_y_pos.resize( nets_size, 0 );
    m_nets_sum_p_inv_x_pos.resize( nets_size, 0 );
    m_nets_sum_p_inv_y_pos.resize( nets_size, 0 );
    m_nets_sum_p_x_neg.resize( nets_size, 0 );
    m_nets_sum_p_y_neg.resize( nets_size, 0 );
    m_nets_sum_p_inv_x_neg.resize( nets_size, 0 );
    m_nets_sum_p_inv_y_neg.resize( nets_size, 0 );

    _cellPotentialNorm.resize( modules_size );

    x        = new double [ 2 * modules_size ];
    _expX    = new double [ 2 * modules_size ];
    _expPins = new double [ 2 * pins_size ];
    x_l      = new double [ 2 * modules_size ];
    x_u      = new double [ 2 * modules_size ];
}
*/
/* Constructor. */
MyNLP::MyNLP( CPlaceDB& db )
    : _potentialGridR( 2 ),
      m_potentialGridSize( -1 ),
      m_targetUtil( 0.9 )
{
    m_lookAheadLegalization = false;
    m_earlyStop = false;
    m_topLevel = false;
    m_lastNLP = false;
    m_useBellPotentialForPreplaced = true;
    
    m_weightWire = 4.0;
    m_incFactor = 2.0;
    m_smoothR = 5;	// Gaussian smooth R
    m_smoothDelta = 1;

    target_nnb = 1.1;	// UNUSE

    //TODO: compute target density  (evenly distribute?)
    
    m_pDB = &db;
    InitModuleNetPinId();
    
    //alpha = 3; //m_pDB->m_rowHeight * 20;
    //_weightDensity = 0.01;
    //_weightWire    = 1000.0;
    //m_potentialGridSize      = 1000;

    // scale between 0 to 10
    const double range = 10.0;
    if( m_pDB->m_coreRgn.right > m_pDB->m_coreRgn.top )
    {
        //printf( "right > top\n" );
        m_posScale = range / m_pDB->m_coreRgn.right;
    }
    else
    {
        //printf( "right < top\n" );
        m_posScale = range / m_pDB->m_coreRgn.top;
    }

    _cellPotentialNorm.resize( m_pDB->m_modules.size() );

    x        = new double [ 2 * m_pDB->m_modules.size() ];
    xBest    = new double [ 2 * m_pDB->m_modules.size() ];
    xBak     = new double [ 2 * m_pDB->m_modules.size() ];	// for stop checking (UNUSE NOW)
    //xBak2    = new double [ 2 * m_pDB->m_modules.size() ];	// for line search  (UNUSE NOW)
    _expX    = new double [ 2 * m_pDB->m_modules.size() ];
    _expPins = new double [ 2 * m_pDB->m_pins.size() ];
    x_l      = new double [ 2 * m_pDB->m_modules.size() ];
    x_u      = new double [ 2 * m_pDB->m_modules.size() ];
    
    m_usePin.resize( m_pDB->m_modules.size() );
    SetUsePin();

    m_nets_sum_exp_xi_over_alpha.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_exp_yi_over_alpha.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_exp_inv_xi_over_alpha.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_exp_inv_yi_over_alpha.resize( m_pDB->m_nets.size(), 0 );

    m_nets_sum_p_x_pos.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_p_y_pos.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_p_inv_x_pos.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_p_inv_y_pos.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_p_x_neg.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_p_y_neg.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_p_inv_x_neg.resize( m_pDB->m_nets.size(), 0 );
    m_nets_sum_p_inv_y_neg.resize( m_pDB->m_nets.size(), 0 );

    grad_wire.resize( 2 * m_pDB->m_modules.size(), 0.0 );
    grad_potential.resize( 2 * m_pDB->m_modules.size(), 0.0 );

    //flpengchange 2015 12.04
    sub_grad_wire.resize(2 * m_pDB->m_modules.size(), 0.0);
    sub_grad_potential.resize(2 * m_pDB->m_modules.size(), 0.0);


    // hdgao added
    part_grad_wire.resize(2 * m_pDB->m_modules.size(), 0.0);
    part_grad_potential.resize(2 * m_pDB->m_modules.size(), 0.0);
    part_grad.resize(2*m_pDB->m_modules.size(),0.0);

    //    for(int i = 0; i < THREAD_size; i++)
    //            save_x_threads[i].resize( 2 * m_pDB->m_modules.size() );


    /*
    m_totalMovableModuleArea = 0;
    m_totalFixedModuleArea = 0;
    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
    if( m_pDB->m_modules[i].m_isFixed == false )
        m_totalMovableModuleArea += m_pDB->m_modules[i].m_area;
    else
    {
        if( m_pDB->m_modules[i].m_isOutCore == false )
        {
        m_totalFixedModuleArea +=
            getOverlap( m_pDB->m_coreRgn.left, m_pDB->m_coreRgn.right,
                m_pDB->m_modules[i].m_x, m_pDB->m_modules[i].m_x + m_pDB->m_modules[i].m_width ) *
            getOverlap( m_pDB->m_coreRgn.bottom, m_pDB->m_coreRgn.top,
                m_pDB->m_modules[i].m_y, m_pDB->m_modules[i].m_y + m_pDB->m_modules[i].m_height );
        }
    }
    }
    if( param.bShow )
    printf( "Total movable area %.0f, fixed area %.0f\n",
        m_totalMovableModuleArea, m_totalFixedModuleArea );
    */
    
}
void MyNLP::MyNLP_copy(MyNLP* source){
    m_lookAheadLegalization = source->m_lookAheadLegalization;
    m_earlyStop = source->m_earlyStop;
    m_topLevel = source->m_topLevel;
    m_lastNLP = source->m_lastNLP;
    m_useBellPotentialForPreplaced = source->m_useBellPotentialForPreplaced;
    m_smoothDelta = source->m_smoothDelta;

    m_weightWire = source->m_weightWire;
    m_incFactor = source->m_incFactor;
    m_smoothR = source->m_smoothR;	// Gaussian smooth R


    target_nnb = source->target_nnb;	// UNUSE

    _potentialGridR = source->_potentialGridR;
    m_potentialGridSize = source->m_potentialGridSize;
    m_targetUtil = source->m_targetUtil;

    m_smoothR = source->m_smoothR;


    m_pDB = new CPlaceDB;
    m_pDB->CPlaceDB_copy( source->m_pDB);

    m_posScale = source->m_posScale;

    _cellPotentialNorm.resize( m_pDB->m_modules.size() );

    x        = new double [ 2 * m_pDB->m_modules.size() ];
    xBest    = new double [ 2 * m_pDB->m_modules.size() ];
    _expX    = new double [ 2 * m_pDB->m_modules.size() ];
    _expPins = new double [ 2 * m_pDB->m_pins.size() ];
    x_l      = new double [ 2 * m_pDB->m_modules.size() ];
    x_u      = new double [ 2 * m_pDB->m_modules.size() ];

    memcpy(x,       source->x,       sizeof(double)*2 * m_pDB->m_modules.size());
    memcpy(xBest,   source->xBest,   sizeof(double)*2 * m_pDB->m_modules.size());
    memcpy(_expX,   source->_expX,   sizeof(double)*2 * m_pDB->m_modules.size());
    memcpy(_expPins,source->_expPins,sizeof(double)*2 * m_pDB->m_pins.size());
    memcpy(x_l,     source->x_l,     sizeof(double)*2 * m_pDB->m_modules.size());
    memcpy(x_u,     source->x_u,     sizeof(double)*2 * m_pDB->m_modules.size());

    m_usePin.resize( m_pDB->m_modules.size() );
    m_usePin.assign(source->m_usePin.begin(),source->m_usePin.end());
    //  SetUsePin();

    _weightWire = source->_weightWire;
    _weightDensity = source->_weightDensity;
    m_ite = source->m_ite;
    m_currentStep = source->m_currentStep;
    m_incFactor = source->m_incFactor;
    m_weightWire = source->m_weightWire;
    _alpha = source->_alpha;
    m_alwaysOverPotential = source->m_alwaysOverPotential;



    m_nets_sum_exp_xi_over_alpha.assign(source->m_nets_sum_exp_xi_over_alpha.begin(),
                                        source->m_nets_sum_exp_xi_over_alpha.end());
    m_nets_sum_exp_yi_over_alpha.assign(source->m_nets_sum_exp_yi_over_alpha.begin(),
                                        source->m_nets_sum_exp_yi_over_alpha.end());
    m_nets_sum_exp_inv_xi_over_alpha.assign(source->m_nets_sum_exp_inv_xi_over_alpha.begin(),
                                            source->m_nets_sum_exp_inv_xi_over_alpha.end());
    m_nets_sum_exp_inv_yi_over_alpha.assign(source->m_nets_sum_exp_inv_yi_over_alpha.begin(),
                                            source->m_nets_sum_exp_inv_yi_over_alpha.end());

    m_nets_sum_p_x_pos.assign(source->m_nets_sum_p_x_pos.begin(),source->m_nets_sum_p_x_pos.end());
    m_nets_sum_p_y_pos.assign(source->m_nets_sum_p_y_pos.begin(),source->m_nets_sum_p_y_pos.end());
 //   m_nets_sum_p_inv_x_pos.assign(source->m_nets_sum_p_inv_x_pos.begin(),source->m_nets_sum_p_inv_x_pos.end());
 //   m_nets_sum_p_inv_y_pos.assign(source->m_nets_sum_p_inv_y_pos.begin(),source->m_nets_sum_p_inv_y_pos.end());
 //   m_nets_sum_p_x_neg.assign(source->m_nets_sum_p_x_neg.begin(),source->m_nets_sum_p_x_neg.end());
 //   m_nets_sum_p_y_neg.assign(source->m_nets_sum_p_y_neg.begin(),source->m_nets_sum_p_y_neg.end());
    m_nets_sum_p_inv_x_neg.assign(source->m_nets_sum_p_inv_x_neg.begin(),source->m_nets_sum_p_inv_x_neg.end());
    m_nets_sum_p_inv_y_neg.assign(source->m_nets_sum_p_inv_y_neg.begin(),source->m_nets_sum_p_inv_y_neg.end());
    _cellPotentialNorm.assign(source->_cellPotentialNorm.begin(),source->_cellPotentialNorm.end());

    m_potentialGridWidth  = source->m_potentialGridWidth;
    m_potentialGridHeight = source->m_potentialGridHeight;
    _potentialRY = source->_potentialRY;
    _potentialRX = source->_potentialRX;

    m_gridPotential.resize(source->m_gridPotential.size());
    m_basePotential.resize(source->m_basePotential.size());
    m_binFreeSpace.resize(source->m_binFreeSpace.size());
    m_basePotentialOri.resize(source->m_basePotentialOri.size());
    m_expBinPotential.resize(source->m_expBinPotential.size());
    m_gridDensity.resize(source->m_gridDensity.size());
    m_gridDensitySpace.resize(source->m_gridDensitySpace.size());
    m_moduleNetPinId.resize(source->m_pDB->m_modules.size());

    for(size_t i = 0; i < m_gridPotential.size();i++)
        m_gridPotential[i].assign(source->m_gridPotential[i].begin(),
                                  source->m_gridPotential[i].end());

    for(size_t i = 0; i < m_basePotential.size();i++)
        m_basePotential[i].assign(source->m_basePotential[i].begin(),
                                  source->m_basePotential[i].end());

    for(size_t i = 0; i < m_binFreeSpace.size();i++)
        m_binFreeSpace[i].assign(source->m_binFreeSpace[i].begin(),
                                 source->m_binFreeSpace[i].end());

    for(size_t i = 0; i < m_basePotentialOri.size();i++)
        m_basePotentialOri[i].assign(source->m_basePotentialOri[i].begin(),
                                     source->m_basePotentialOri[i].end());

    for(size_t i = 0; i < m_expBinPotential.size();i++)
        m_expBinPotential[i].assign(source->m_expBinPotential[i].begin(),
                                    source->m_expBinPotential[i].end());

    for(size_t i = 0; i < m_gridDensity.size();i++)
        m_gridDensity[i].assign(source->m_gridDensity[i].begin(),
                                source->m_gridDensity[i].end());

    for(size_t i = 0; i < m_gridDensitySpace.size();i++)
        m_gridDensitySpace[i].assign(source->m_gridDensitySpace[i].begin(),
                                     source->m_gridDensitySpace[i].end());
    for(size_t i = 0; i < m_moduleNetPinId.size();i++)
        m_moduleNetPinId[i].assign(source->m_moduleNetPinId[i].begin(),
                                   source->m_moduleNetPinId[i].end());

    grad_wire.resize( 2 * m_pDB->m_modules.size(), 0.0 );
    grad_potential.resize( 2 * m_pDB->m_modules.size(), 0.0 );

    //flpengchange 2015 12.04
    sub_grad_wire.resize(2 * m_pDB->m_modules.size(), 0.0);
    sub_grad_potential.resize(2 * m_pDB->m_modules.size(), 0.0);


    // hdgao added
    part_grad_wire.resize(2 * m_pDB->m_modules.size(), 0.0);
    part_grad_potential.resize(2 * m_pDB->m_modules.size(), 0.0);
    part_grad.resize(2*m_pDB->m_modules.size(),0.0);


    m_alwaysOverArea    = source->m_alwaysOverArea;
    m_gridDensityWidth  = source->m_gridDensityWidth;
    m_gridDensityHeight = source->m_gridDensityHeight;
    m_gridDensityTarget = source->m_gridDensityTarget;

    PL_step = source->PL_step;








}

MyNLP::~MyNLP()
{
    delete [] x;
    delete [] xBest;
    //delete [] xBak;
    //delete [] xBak2;
    delete [] _expX;
    delete [] _expPins;
    delete [] x_l;
    delete [] x_u;
}

void MyNLP::SetUsePin()
{
    //printf( "row height = %f\n", m_pDB->m_rowHeight );
    int effectivePinCount = 0;
    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
        bool usePin = false;
        for( unsigned int p=0; p<m_pDB->m_modules[i].m_pinsId.size(); p++ )
        {
            int pinId = m_pDB->m_modules[i].m_pinsId[p];

            if( m_pDB->m_pins[pinId].xOff != 0.0 || m_pDB->m_pins[pinId].yOff != 0.0 )
            {
                usePin = true;
                break;
            }
        }
        if( usePin )
            effectivePinCount++;
        m_usePin[i] = usePin;
    }
    //   if( param.bShow )
    //	printf( "Effective Pin # = %d\n", effectivePinCount );
}



bool MyNLP::MySolve(double wWire,
                    double target_density,
                    int currentLevel,	// for plotting
                    bool noRelaxSmooth,
                    threadpool_t **pool)
{

    double time_start = seconds();
    //printf( "Start Optimization \n" );
    assert( _potentialGridR > 0 );
    double avgCellSize = 0;
    int count = 0;
    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
        if( m_pDB->m_modules[i].m_isFixed == false )
        {
            avgCellSize += sqrt( m_pDB->m_modules[i].m_area );
            count++;
        }
    }
    avgCellSize /= count;
    //int maxGridSize = static_cast<int>( (m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left) / avgCellSize * 0.95 );
    int maxGridSize = static_cast<int>( (m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left) / avgCellSize );
    
    
    if( m_potentialGridSize <= 0 )
        //m_potentialGridSize = static_cast<int>( sqrt( m_pDB->m_modules.size() ) * 0.8 );
        //    m_potentialGridSize = static_cast<int>( sqrt(static_cast<double>( m_pDB->m_modules.size()) ) * 0.95 );
        m_potentialGridSize = static_cast<int>( sqrt(static_cast<double>( m_pDB->m_modules.size()) ) * 0.8 );//original 0.8




    if( param.bShow )
        printf( "step scale = %f\n", param.step );
    
    int n, m, nnz_jac_g, nnz_h_lag;
    get_nlp_info( n, m, nnz_jac_g, nnz_h_lag );
    get_bounds_info( n, x_l, x_u, m, nullptr, nullptr );
    get_starting_point( n, true, x, false, nullptr, nullptr, m, false, nullptr );
    BoundX( n, x, x_l, x_u );

    m_ite = 0;
    bool isLegal = false;
    assert( param.dLpNorm_P > 0 );
    while( true )
    {
        //_alpha = 0.5 * m_potentialGridWidth; // according to APlace ispd04
        //_alpha = ( m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left ) * 0.005;	// as small as possible
        _alpha = param.dLpNorm_P;
        //printf( "GRID = %d (alpha = %f) width = %.2f\n", m_potentialGridSize, _alpha, ( m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left )/m_potentialGridSize );
        if( param.bShow )
            printf( "GRID = %d  (width = %.2f)\n", m_potentialGridSize, ( m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left )/m_potentialGridSize );

        if( m_topLevel )
            m_lastNLP = true;
        else
            m_lastNLP = false;
        isLegal = ParallelSolve(wWire, target_density, currentLevel,pool );
       // isLegal = PLSolve_smooth( wWire, target_density, currentLevel,pool );
         //      isLegal = SLSolve_smooth( wWire, target_density, currentLevel );
         //    isLegal = GoSolve( wWire, target_density, currentLevel );
        break;

        if( !m_topLevel )
            break;

        m_potentialGridSize *= 2;	    // finner the grid when top level

        if( m_potentialGridSize > maxGridSize )
            break;
    }

    if( param.bShow )
    {
        printf( "HPWL = %.0f\n", m_pDB->CalcHPWL() );
        printf( "\nLevel Time        = %.2f sec = %.2f min\n",
                double(seconds()-time_start), double(seconds()-time_start)/60.0  );
        printf( "Time Sum-Exp      = %.1f sec\n", time_sum_exp );
        printf( "Time up potential = %.1f sec\n", time_up_potential );
        printf( "Time eval_f       = %.2f sec = (WL) %.2f + (P)%.2f\n", time_f, time_wl, time_update_grid );
        printf( "Time eval_grad_f  = %.1f sec = (gradWL) %.1f + (gradP) %.1f\n",
                time_grad_f, time_grad_wl, time_grad_potential );
    }

    return isLegal;
}


bool MyNLP::GoSolve( double wWire, 
                     double target_density,
                     int currentLevel	// for plotting
                     )
{
    double givenTargetUtil = m_targetUtil; // for look ahead legalization
    m_currentStep = param.step;
    
    m_targetUtil += 0.05;
    if( m_targetUtil > 1.0 )
        m_targetUtil = 1.0;
    
    
    double time_start = seconds();
    char filename[100];	    // for gnuplot
    

    size_t n = 2 * m_pDB->m_modules.size();
    
    double designUtil = m_pDB->m_totalMovableModuleArea / m_pDB->m_totalFreeSpace;
    if( param.bShow )
        printf( "INFO: Design utilization: %f\n", designUtil );
    if( m_targetUtil > 0 )  // has user-defined target utilization
    {
        // This part is very important for ISPD-06 Placement Contest.
        double lowest = designUtil + 0.05;
        if( m_targetUtil < lowest )
        {
            if( param.bShow )
            {
                printf( "WARNING: Target utilization (%f) is too low\n", m_targetUtil );
                printf( "         Set target utilization to %f\n", lowest );
            }
            m_targetUtil = lowest;
        }
    }
    else // no given utilization
    {
        printf( "No given target utilization. Distribute blocks evenly.\n" );
        m_targetUtil = designUtil + 0.05;
        if( m_targetUtil > 1.0 )
            m_targetUtil = 1.0;
    }
    if( param.bShow )
        printf( "DBIN: Target utilization: %f\n", m_targetUtil );
    
    double* grad_f = new double [ 2 * m_pDB->m_modules.size() ];
    double* last_grad_f = new double [ 2 * m_pDB->m_modules.size() ]; // for computing CG-direction;
    memset( grad_f, 0, sizeof(double)*2*m_pDB->m_modules.size() );
    memset( last_grad_f, 0, sizeof(double)*2*m_pDB->m_modules.size() );

    /*for( unsigned int i=0; i<2*m_pDB->m_modules.size(); i++ )
    {
    grad_f[i] = 0.0;
    last_grad_f[i] = 0.0;
    }*/

    CreatePotentialGrid();   // create potential grid according to "m_potentialGridSize"
    int densityGridSize = 10;	// 1% chip area
    //int densityGridSize = m_potentialGridSize / 3;	    // not good in big3
    CreateDensityGrid( densityGridSize );	// real density: use 1% area
    UpdateDensityGridSpace( n, x );
    UpdatePotentialGridBase( x );		// init exp potential for each bin, also update ExpBin
    
#if 1
    // gaussian smoothing for base potential
    GaussianSmooth smooth;
    int r = m_smoothR;
    smooth.Gaussian2D( r, 6*r+1 );
    smooth.Smooth( m_basePotential );
    m_basePotentialOri = m_basePotential;
    sprintf( filename, "gbase%d.dat", currentLevel );
    OutputPotentialGrid( filename );
#endif 


    // TEST
    if( m_smoothDelta == 1 )
    {
        if( param.bShow )
        {
            sprintf( filename, "gbase%d-more.dat", currentLevel );
            printf( "generate %s...\n", filename );
            fflush( stdout );
        }

        vector< vector< double > > moreSmooth = m_basePotential;
        r = m_smoothR * 6;
        int kernel_size = 5*r;
        if( kernel_size % 2 == 0 )
            kernel_size++;
        smooth.Gaussian2D( r, kernel_size );
        smooth.Smooth( moreSmooth );

        if( param.bShow )
        {
            swap( moreSmooth, m_basePotential );
            OutputPotentialGrid( filename );
            swap( moreSmooth, m_basePotential );
        }

        // merge base and moreSmooth
        double binArea = m_potentialGridWidth * m_potentialGridHeight;
        double halfBinArea = binArea / 2;
        int changeCount = 0;
        for( unsigned int i=0; i<moreSmooth.size(); i++ )
        {
            for( unsigned int j=0; j<moreSmooth[i].size(); j++ )
            {
                double free = binArea - m_basePotential[i][j];
                if( free < 1e-4 )	// no space
                {
                    if( moreSmooth[i][j] > halfBinArea )
                    {
                        m_basePotential[i][j] += moreSmooth[i][j] - halfBinArea;
                        changeCount++;
                    }
                }
            }
        }

        if( param.bShow )
        {
            printf( "change %d\n", changeCount );
            sprintf( filename, "gbase%d-more-merge.dat", currentLevel );
            OutputPotentialGrid( filename );
        }
    }


    if( m_smoothDelta > 1.0 )
        SmoothPotentialBase( double(m_smoothDelta) );   // also update ExpBin
    
    UpdateExpBinPotential( m_targetUtil );

#if 1 
    if( param.bShow )
    {
        sprintf( filename, "base%d.dat", currentLevel );
        OutputPotentialGrid( filename );
        // TEST
        /*for( int delta=1; delta<=10; delta++ )
      {
      SmoothPotentialBase( (double)delta );
      sprintf( filename, "base%d-%d.dat", currentLevel, delta );
      OutputPotentialGrid( filename );
      }*/
    }
#endif


    assert( m_targetUtil > 0 );
    
    // wirelength
    UpdateExpValueForEachCell( n, x, _expX, _alpha );
    UpdateExpValueForEachPin( n, x, _expPins, _alpha );
    UpdateNetsSumExp( x, _expX );
    totalWL = GetWL( n, x, _expX, _alpha );

    // density
    UpdatePotentialGrid( x );
    UpdateDensityGrid( n, x );
    density = GetDensityPanelty();

    // 2006-02-22 weight (APlace ICCAD05)
    _weightWire = 1.0;
    eval_grad_f( n, x, _expX, true, grad_f );
    double totalWireGradient = 0;
    double totalPotentialGradient = 0;

    // TODO: truncation?
    AdjustForce( n, x, grad_wire, grad_potential );

    for( int i=0; i<n; i++ )
    {
        totalWireGradient      += fabs( grad_wire[i] );
        totalPotentialGradient += fabs( grad_potential[i] );
        //totalWireGradient      += grad_wire[i] * grad_wire[i];
        //totalPotentialGradient += grad_potential[i] * grad_potential[i];
    }
    //totalWireGradient = sqrt( totalWireGradient );
    //totalPotentialGradient = sqrt( totalPotentialGradient );

    /*
    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
    totalWireGradient += sqrt( grad_wire[2*i] * grad_wire[2*i] + grad_wire[2*i+1] * grad_wire[2*i+1] );
    totalPotentialGradient += sqrt( grad_potential[2*i] * grad_potential[2*i] + grad_potential[2*i+1] * grad_potential[2*i+1] );
    }
    */
    
    //_weightDensity = 1.0 / ( 0.5 * totalPotentialGradient / totalWireGradient );  // APlace ICCAD-05
    //_weightDensity = 1.0 / ( wWire * totalPotentialGradient / totalWireGradient );
    
    _weightDensity = 1.0;
    _weightWire = wWire * totalPotentialGradient / totalWireGradient;
    
    //    printf( " INIT: LogSumExp WL= %.0f, gradWL= %.0f\n", totalWL, totalWireGradient );
    //    printf( " INIT: DensityPenalty= %.0f, gradPenalty= %.0f\n", density, totalPotentialGradient );
    
    int maxIte = 50;	// max outerIte

    //    printf( "       DenGridSize   = %d\n", densityGridSize );
    //    printf( "       GridSize      = %d * %d (w %.3f h %.3f)\n", m_potentialGridSize, m_potentialGridSize, m_potentialGridWidth, m_potentialGridHeight );
    //    printf( "       GridLen       = w %.3f h %.3f\n", m_potentialGridWidth, m_potentialGridHeight );
    //    printf( "       alpha         = %f\n", _alpha );
    //    printf( "       maxOuterIter  = %d\n", maxIte );
    //    printf( "       weightDensity = %g\n", _weightDensity );
    //    printf( "       weightWL      = %g\n", _weightWire );
    //    printf( "       targetNNB     = %f\n", target_nnb );
    //    printf( "       stopStepSize  = %f\n", stopStepSize );
    

    
    bool newDir = true;
    double obj_value;
    double beta;	// determined by CG
    eval_f( n, x, _expX, true, obj_value );
    eval_grad_f( n, x, _expX, true, grad_f );

    double nnb_real = GetNonZeroDensityGridPercent();
    UpdateDensityGrid( n, x );
    double maxDen = GetMaxDensity();
    double totalOverDen = GetTotalOverDensity();
    double totalOverDenLB = GetTotalOverDensityLB();
    double totalOverPotential = GetTotalOverPotential();
    
    //printf( "INIT f = %f\n", obj_value );

    if( param.bShow )
    {
        printf( " %d-%2d HPWL= %.0f\tDen= %.2f %.2f %.2f %.2f NNB= %.2f Dcost= %4.1f%%  WireW= %.0f ",
                currentLevel, m_ite, m_pDB->CalcHPWL(),
                maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                nnb_real,
                density * _weightDensity / obj_value * 100.0, _weightWire
                );
    }
    else
    {
        printf( " %d-%2d HPWL= %.0f \t",
                currentLevel, m_ite, m_pDB->CalcHPWL()
                );
    }
    fflush( stdout );
    if( param.bShow )
    {
        sprintf( filename, "fig%d-%d.plt", currentLevel, m_ite );
        m_pDB->OutputGnuplotFigure( filename, false, false );
    }

    double lastTotalOver = 0;
    double lastTotalOverPotential = DBL_MAX;
    double over = totalOverDen;
    int totalIte = 0;

    bool hasBestLegalSol = false;
    double bestLegalWL = DBL_MAX;
    int lookAheadLegalCount = 0;
    double totalLegalTime = 0.0;
    //double norm_move = 0.0;
    
    bool startDecreasing = false;
    
    int checkStep = 5;
    int outStep = 25;
    if( param.bShow == false )
        outStep = INT_MAX;

    int LALnoGoodCount = 0;
    //  now optimization begin
    for( int ite=0; ite<maxIte; ite++ )
    {
        m_ite++;
        int innerIte = 0;
        double old_obj = DBL_MAX;
        double last_obj_value = DBL_MAX;

        m_currentStep = param.step;

        newDir = true;
        while( true )	// inner loop, minimize "f"
        {
            innerIte++;
            swap( last_grad_f, grad_f );    // save for computing the congujate gradient direction
            eval_grad_f( n, x, _expX, true, grad_f );
            AdjustForce( n, x, grad_f );

            if( innerIte % checkStep == 0 )
            {
                old_obj = last_obj_value;    // backup the old value
                eval_f( n, x, _expX, true, obj_value );
                last_obj_value = obj_value;
            }

#if 1
            // Output solving progress
            if( innerIte % outStep == 0 && innerIte != 0 )
            {
                if( innerIte % checkStep != 0 )
                    eval_f( n, x, _expX, true, obj_value );
                printf( "\n\t  (%4d): f= %.10g\tstep= %.6f \t %.1fm ",
                        innerIte,
                        obj_value,
                        stepSize,
                        double(seconds()-time_start)/60.0
                        );
                /*printf( "\n\t  (%4d): f= %.8g\tstep= %.6f %.2f \t %.1fm ",
            innerIte,
            obj_value,
            stepSize,
            sqrt(norm_move/n),
            double(seconds()-time_start)/60.0
              );*/
                /*UpdateBlockPosition( x );   // update to placeDB
          UpdateDensityGrid( n, x );
          double nnb_real = GetNonZeroDensityGridPercent();
          maxDen = GetMaxDensity();
          totalOverDen = GetTotalOverDensity();
          printf( "\n %4d HPWL= %0.4g\tDen= %.3f %.3f NNB= %.3f LTime= %.1fm ",
          innerIte, m_pDB->CalcHPWL(),
          maxDen, totalOverDen,
          nnb_real,
          double(seconds()-time_start)/60.0  );*/
                fflush( stdout );
            }
#endif

            if( innerIte % checkStep == 0 )
            {
                printf( "." );
                fflush( stdout );

                if( innerIte % 2 * checkStep == 0 )
                {
                    UpdateBlockPosition( x );   // update to placeDB
                    if( m_pDB->CalcHPWL() > bestLegalWL )   // gWL > LAL-WL
                    {
                        printf( "X\n" );
                        fflush( stdout );
                        break;
                    }
                }

                UpdateDensityGrid( n, x );  // find the exact bin density
                totalOverDen = GetTotalOverDensity();
                totalOverDenLB = GetTotalOverDensityLB();
                totalOverPotential = GetTotalOverPotential();

                lastTotalOver = over;
                over = min( totalOverPotential, totalOverDen ); // TEST

                if( !startDecreasing
                        && over < lastTotalOver
                        && ite >= 1
                        && innerIte >= 6 )
                {
                    printf( ">>" );
                    fflush( stdout );
                    startDecreasing = true;
                }

                // 2005-03-11: meet the constraint
                if( startDecreasing && over < target_density )
                    break;

                // Cannot further improve
                if( obj_value >= param.precision * old_obj )
                {
                    break;
#if 0
                    if( m_currentStep < 0.2 )
                        break;
                    else
                    {
                        //m_currentStep *= 0.618;
                        m_currentStep *= 0.6666;
                        printf( "*" );
                        newDir = true;
                        printf( "\n\t  (%4d): f= %.10g\tstep= %.6f \t %.1fm ",
                                innerIte,
                                obj_value,
                                stepSize,
                                double(seconds()-time_start)/60.0
                                );
                    }
#endif
                }
            }


            // Calculate d_k (conjugate gradient method)
            if( newDir == true )
            {
                // gradient direction
                newDir = false;
                for( int i=0; i<n; i++ )
                    grad_f[i] = -grad_f[i];
            }
            else
            {
                // conjugate gradient direction
                FindBeta( n, grad_f, last_grad_f, beta );
                for( int i=0; i<n; i++ )
                    grad_f[i] = -grad_f[i] + beta * last_grad_f[i];
            }


            // Calculate a_k (step size)
            LineSearch( n, x, grad_f, stepSize );

            // Update X. (x_{k+1} = x_{k} + \alpha_k * d_k)
            double move;
            for( int i=0; i<n; i++ )
            {
                move = grad_f[i] * stepSize;
                x[i] += move;
            }

            /*
        norm_move = 0.0;
        for( int i=0; i<n; i++ )
        {
        double move = grad_f[i] * stepSize;
        norm_move += move * move;
        }
        */

            BoundX( n, x, x_l, x_u );
            double time_used = seconds();
            UpdateExpValueForEachCell( n, x, _expX, _alpha );
            UpdateExpValueForEachPin( n, x, _expPins, _alpha );
            UpdateNetsSumExp( x, _expX );
            time_grad_wl += seconds() - time_used;
            UpdatePotentialGrid( x );

        }// while(true) end   inner loop

        if( param.bShow )
        {
            printf( "%d\n", innerIte );
            fflush( stdout );
        }
        else
            printf( "\n" );
        totalIte += innerIte;

        UpdateDensityGrid( n, x );
        double nnb_real = GetNonZeroDensityGridPercent();
        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();

        UpdateBlockPosition( x );   // update to placeDB

#if 1
        if( param.bShow )
        {
            // output figures
            sprintf( filename, "fig%d-%d.plt", currentLevel, m_ite );
            m_pDB->OutputGnuplotFigure( filename, false, false );	// it has "CalcHPWL()"
            if( m_topLevel ) // debugging
            {
                sprintf( filename, "fig%d-%d.pl", currentLevel, m_ite );
                m_pDB->OutputPL( filename );
            }
            //if( m_potentialGridSize < 30 )  PrintPotentialGrid();

            sprintf( filename, "grid%d-%d.dat", currentLevel, m_ite );
            OutputPotentialGrid( filename );
            sprintf( filename, "den%d-%d.dat", currentLevel, m_ite );
            OutputDensityGrid( filename );

            sprintf( filename, "util%d-%d.dat", currentLevel, m_ite );
            CPlaceBin placeBin( *m_pDB );
            placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
            placeBin.OutputBinUtil( filename );
        }
#endif

        if( param.bShow )
        {
            printf( " %d-%2d HPWL= %.0f\tDen= %.2f %.4f %.4f %.4f NNB= %.2f LTime= %.1fm Dcost= %4.1f%% WireW= %.0f ",
                    currentLevel, m_ite, m_pDB->CalcHPWL(),
                    maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                    nnb_real,
                    double(seconds()-time_start)/60.0,
                    density*_weightDensity /obj_value * 100.0,
                    _weightWire
                    );
        }
        else
        {
            printf( " %d-%2d HPWL= %.f\tLTime= %.1fm ",
                    currentLevel, m_ite, m_pDB->CalcHPWL(),
                    double(seconds()-time_start)/60.0
                    );
        }
        fflush( stdout );


#if 1
        // 2006-03-06 (CAUTION! Do not use look-ahead legalization when dummy block exists.
        // TODO: check if there is dummy block (m_modules[].m_isDummy)
        if( m_ite >= 2 && m_lookAheadLegalization && over < target_density+0.10 )
            //if( startDecreasing && m_lookAheadLegalization ) // test
        {
            UpdateBlockPosition( x );   // update to placeDB
            double hpwl = m_pDB->CalcHPWL();
            if( hpwl > bestLegalWL )
            {
                printf( "Stop. Good enough.1\n" );
                break;
            }

            lookAheadLegalCount++;
            double oldWL = hpwl;
            CTetrisLegal legal(*m_pDB);

            //bool bMacroShifter = legal.MacroShifter( 10, false );
            //if( false == bMacroShifter )
            //	printf( "MACRO SHIFTER FAILED!\n" );

            double scale = 0.85;
            if( givenTargetUtil < 1.0 && givenTargetUtil > 0 )
                scale = 0.9;

            double legalStart = seconds();
            bool bLegal = legal.Solve( givenTargetUtil, false, false, scale );
            double legalTime = seconds() - legalStart;
            totalLegalTime += legalTime;
            if( param.bShow )
                printf( "LAL Time: %.2f\n", legalTime );
            if( bLegal )
            {
                double WL = m_pDB->GetHPWLdensity( givenTargetUtil );
                if( param.bShow )
                    m_pDB->ShowDensityInfo();
                if( WL < bestLegalWL )
                {
                    // record the best legal solution
                    LALnoGoodCount = 0;
                    if( param.bShow )
                        printf( "SAVE BEST! (HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n",
                                m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    bestLegalWL = WL;
                    hasBestLegalSol = true;
                    for( int i=0; i<(int)m_pDB->m_modules.size(); i++ )
                    {
                        xBest[2*i] = m_pDB->m_modules[i].m_cx;
                        xBest[2*i+1] = m_pDB->m_modules[i].m_cy;
                    }
                }
                else
                {
                    if( param.bShow )
                        printf( "(HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n", m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    if( (WL-oldWL)/oldWL < 0.075 )
                    {
                        if( param.bShow )
                            printf( "Stop. Good enough 2\n" );
                        break;
                    }
                    LALnoGoodCount++;
                    if( LALnoGoodCount >= 2 )
                        break;
                }
            }
        }
#endif	

        if( ite >= 2 )
        {
            if( startDecreasing && over < target_density )
            {
                printf( "Meet constraint!\n" );
                break;
            }

            // cannot reduce totalOverPotential
            if( ite > 3 && totalOverPotential > lastTotalOverPotential &&
                    totalOverPotential < 1.4 )
            {
                printf( "Cannot further reduce!\n" );
                break;
            }
        }

        _weightWire /= m_incFactor;
        lastTotalOverPotential = totalOverPotential;

    }// outer loop


    // 2006-03-06 (donnie)
    if( hasBestLegalSol )
        memcpy( x, xBest, sizeof(double)*n );
    UpdateBlockPosition( x );

    if( lookAheadLegalCount > 0 && param.bShow )
    {
        printf( "LAL: Total Count: %d\n", lookAheadLegalCount );
        printf( "LAL: Total CPU: %.2f\n", totalLegalTime );
        sprintf( filename, "util-global.dat" );
        CPlaceBin placeBin( *m_pDB );
        placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
        placeBin.OutputBinUtil( filename );
    }
    
    static int allTotalIte = 0;
    allTotalIte += totalIte;

    if( param.bShow )
    {
        m_pDB->ShowDensityInfo();
        printf( "\nLevel Ite %d   Total Ite %d\n", totalIte, allTotalIte );
    }
    
    delete [] grad_f;
    delete [] last_grad_f;	// for computing conjugate gradient direction
    
    return hasBestLegalSol;
}


void MyNLP::FindBeta( const int& n, const double* grad_f, const double* last_grad_f, double& beta )
{
    // Polak-Ribiere foumula from APlace journal paper
    // NOTE:
    //   g_{k-1} = -last_grad_f
    //   g_k     = grad_f

    double l2norm = 0;
    for( int i=0; i<n; i++ )
        l2norm += last_grad_f[i] * last_grad_f[i];

    double product = 0;
    for( int i=0; i<n; i++ )
        product += grad_f[i] * ( grad_f[i] + last_grad_f[i] );	// g_k^T ( g_k - g_{k-1} )
    beta = product / l2norm;
    param.sum_step += beta;
    param.count ++;
}


void MyNLP::BoundX( const int& n, double* x, double* x_l, double* x_h, const int& i )
{
    if( x[i] < x_l[i] )      x[i] = x_l[i];
    else if( x[i] > x_h[i] )	x[i] = x_h[i];
}


void MyNLP::BoundX( const int& n, double* x, double* x_l, double* x_h )
{
    for( int i=0; i<n; i++ )
    {
        if( x[i] < x_l[i] )             x[i] = x_l[i];
        else if( x[i] > x_h[i] )	x[i] = x_h[i];
    }
}

void MyNLP::AdjustForce( const int& n, const double* x, double* grad_f )
{

    double totalGrad = 0;
    int size = n/2;
    for( int i=0; i<size; i++ )
    {
        double value = grad_f[2*i] * grad_f[2*i] + grad_f[2*i+1] * grad_f[2*i+1];
        totalGrad += value;
    }
    double avgGrad = sqrt( totalGrad / size );

    // Do truncation
    double expMaxGrad = avgGrad * truncationFactor;	// x + y
    double expMaxGradSquare = expMaxGrad * expMaxGrad;
    for( int i=0; i<size; i++ )
    {
        double valueSquare = ( grad_f[2*i]*grad_f[2*i] + grad_f[2*i+1]*grad_f[2*i+1] );
        if( valueSquare > expMaxGradSquare )
        {
            double value = sqrt( valueSquare );
            grad_f[2*i]   = grad_f[2*i]   * expMaxGrad / value;
            grad_f[2*i+1] = grad_f[2*i+1] * expMaxGrad / value;
        }
    }
}
void MyNLP::AdjustForce( const int& n, const double* x, vector<double>& grad_f )
{

    double totalGrad = 0;
    int size = n/2;
    for( int i=0; i<size; i++ )
    {
        double value = grad_f[2*i] * grad_f[2*i] + grad_f[2*i+1] * grad_f[2*i+1];
        totalGrad += value;
    }
    double avgGrad = sqrt( totalGrad / size );

    // Do truncation
    double expMaxGrad = avgGrad * truncationFactor;	// x + y
    double expMaxGradSquare = expMaxGrad * expMaxGrad;
    for( int i=0; i<size; i++ )
    {
        double valueSquare = ( grad_f[2*i]*grad_f[2*i] + grad_f[2*i+1]*grad_f[2*i+1] );
        if( valueSquare > expMaxGradSquare )
        {
            double value = sqrt( valueSquare );
            grad_f[2*i]   = grad_f[2*i]   * expMaxGrad / value;
            grad_f[2*i+1] = grad_f[2*i+1] * expMaxGrad / value;
        }
    }
}
void MyNLP::AdjustForce_2D( const int& n, const double* x, vector<double>& grad_f )
{

    double totalGrad1 = 0,totalGrad2 = 0;
    int size = n/2;
    for( int i=0; i<size; i++ )
    {
        double value1 = grad_f[2*i]   * grad_f[2*i];
        double value2 = grad_f[2*i+1] * grad_f[2*i+1];
        totalGrad1 += value1;
        totalGrad2 += value2;
    }

    double avgGrad1 = sqrt( totalGrad1 / size );
    double avgGrad2 = sqrt( totalGrad2 / size );

    // Do truncation
    double expMaxGrad1 = avgGrad1 * truncationFactor;	// x + y
    double expMaxGrad2 = avgGrad2 * truncationFactor;
    double expMaxGradSquare1 = expMaxGrad1 * expMaxGrad1;
    double expMaxGradSquare2 = expMaxGrad2 * expMaxGrad2;
    for( int i=0; i<size; i++ )
    {
        double valueSquare1 = ( grad_f[2*i]   * grad_f[2*i]);
        double valueSquare2 = ( grad_f[2*i+1] * grad_f[2*i+1] );
        if( valueSquare1 > expMaxGradSquare1 )
        {
            double value1 = sqrt( valueSquare1 );
            grad_f[2*i]   = grad_f[2*i]   * expMaxGrad1 / value1;
        }
        if( valueSquare2 > expMaxGradSquare2 )
        {
            double value2 = sqrt( valueSquare2 );
            grad_f[2*i+1] = grad_f[2*i+1] * expMaxGrad2 / value2;
        }
    }
}


void MyNLP::AdjustForce( const int& n, const double* x, vector<double> grad_wl, vector<double> grad_potential )
{
    double totalGrad = 0;
    int size = n/2;
    for( int i=0; i<size; i++ )
    {
        double value =
                (grad_wl[2*i] + grad_potential[2*i]) * (grad_wl[2*i] + grad_potential[2*i]) +
                (grad_wl[2*i+1] + grad_potential[2*i+1]) * (grad_wl[2*i+1] + grad_potential[2*i+1]);
        totalGrad += value;
    }
    double avgGrad = sqrt( totalGrad / size );

    // Do truncation
    double expMaxGrad = avgGrad * truncationFactor;	// x + y
    double expMaxGradSquare = expMaxGrad * expMaxGrad;
    for( int i=0; i<size; i++ )
    {
        double valueSquare =
                (grad_wl[2*i] + grad_potential[2*i]) * (grad_wl[2*i] + grad_potential[2*i]) +
                (grad_wl[2*i+1] + grad_potential[2*i+1]) * (grad_wl[2*i+1] + grad_potential[2*i+1]);
        if( valueSquare > expMaxGradSquare )
        {
            double value = sqrt( valueSquare );
            grad_wl[2*i]   = grad_wl[2*i]   * expMaxGrad / value;
            grad_wl[2*i+1] = grad_wl[2*i+1] * expMaxGrad / value;
            grad_potential[2*i]   = grad_potential[2*i]   * expMaxGrad / value;
            grad_potential[2*i+1] = grad_potential[2*i+1] * expMaxGrad / value;
        }
    }
}


void MyNLP::LineSearch( const int& n, /*const*/ double* x, double* grad_f, double& stepSize )
{
    int size = n / 2;
    double totalGrad = 0;
    for( int i=0; i<n; i++ )
        totalGrad += grad_f[i] * grad_f[i];
    double avgGrad = sqrt( totalGrad / size );
    stepSize = m_potentialGridWidth / avgGrad * m_currentStep;
    
    return;
}
void MyNLP::LineSearch( const int& n, /*const*/ double* x, vector<double>& grad_f, double& stepSize )
{
    int size = n / 2;
    double totalGrad = 0;
    for( int i=0; i<n; i++ )
        totalGrad += grad_f[i] * grad_f[i];
    double avgGrad = sqrt( totalGrad / size );
    stepSize = m_potentialGridWidth / avgGrad * m_currentStep;
    param.stepRecord.push_back(stepSize);
    param.sum_step += stepSize;
    param.count++;
    return;
}

bool MyNLP::get_nlp_info(int& n, int& m, int& nnz_jac_g, 
                         int& nnz_h_lag/*, IndexStyleEnum& index_style*/)
{
    //printf( "*** get_nlp_info() ***\n" );

    n = m_pDB->m_modules.size() * 2;
    
    //printf( "alpha = %f, wireWeigtht = %f, densityWeight = %f\n",
    //    alpha, _weightWire, _weightDensity );
    
    m = 0;	    // no constraint
    nnz_jac_g = 0;  // 0 nonzeros in the jacobian since no constraint

    /*
    // calculate nnz_h
    set< pair<int,int> > cell_pair;
    for( unsigned int i=0; i<m_pDB->m_nets.size(); i++ )
    {
    if( m_pDB->m_nets[i].size() <= 1 )
        continue;

    for( unsigned int first=0; first<m_pDB->m_nets[i].size()-1; first++ )
    {
        for( unsigned int second=1; second<m_pDB->m_nets[i].size(); second++ )
        {
        if( m_pDB->m_nets[i][first] == m_pDB->m_nets[i][second] )
            continue;

        assert( m_pDB->m_nets[i][first] < (int)m_pDB->m_pins.size() );
        assert( m_pDB->m_nets[i][second] < (int)m_pDB->m_pins.size() );

        int blockFirst = m_pDB->m_pins[ m_pDB->m_nets[i][first] ].moduleId;
        int blockSecond = m_pDB->m_pins[ m_pDB->m_nets[i][second] ].moduleId;

        if( blockFirst >= blockSecond )
            cell_pair.insert( pair<int,int>( blockFirst, blockSecond ) );
        else
            cell_pair.insert( pair<int,int>( blockSecond, blockFirst ) );
        }
    }
    }
    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
    if( m_pDB->m_modules[i].m_isFixed == false )
    {
        // diagonal terms
        cell_pair.insert( pair<int,int>( i, i ) );
    }
    }
    _cellPair.reserve( cell_pair.size() );
    copy( cell_pair.begin(), cell_pair.end(), back_inserter( _cellPair ) );
    nnz_h_lag = _cellPair.size() * 2;	    // dx1dx2 & dy1dy2
    */
    //printf( "    nnz_h_lag = %d\n", nnz_h_lag );

    // We use the standard fortran index style for row/col entries
    //index_style = TNLP::C_STYLE;

    return true;
}

bool MyNLP::get_bounds_info(int n, double* x_l, double* x_u,
                            int m, double* g_l, double* g_u)
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    assert(n == (int)m_pDB->m_modules.size() * 2);
    assert(m == 0);

    //printf( "get_bounds_info\n" );

    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
        if( m_pDB->m_modules[i].m_isFixed )
        {
            x_l[2*i] = m_pDB->m_modules[i].m_cx;
            x_u[2*i] = m_pDB->m_modules[i].m_cx;
            x_l[2*i+1] = m_pDB->m_modules[i].m_cy;
            x_u[2*i+1] = m_pDB->m_modules[i].m_cy;
        }
        else
        {
            x_l[2*i]   = m_pDB->m_coreRgn.left   + m_pDB->m_modules[i].m_width  * 0.5;
            x_u[2*i]   = m_pDB->m_coreRgn.right  - m_pDB->m_modules[i].m_width  * 0.5;
            x_l[2*i+1] = m_pDB->m_coreRgn.bottom + m_pDB->m_modules[i].m_height * 0.5;
            x_u[2*i+1] = m_pDB->m_coreRgn.top    - m_pDB->m_modules[i].m_height * 0.5;
        }
    }

    return true;
}

bool MyNLP::get_starting_point(int n, bool init_x, double* x,
                               bool init_z, double* z_L, double* z_U,
                               int m, bool init_lambda,
                               double* lambda)
{
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    //printf( "get_starting_point\n" );

    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
        x[2*i]   = m_pDB->m_modules[i].m_cx;
        x[2*i+1] = m_pDB->m_modules[i].m_cy;
    }

    return true;
}

void MyNLP::UpdateExpValueForEachCell( const int& n, const double* x, double* expX, const double& inAlpha )
{
    for( int i=0; i<n; i++ )
    {
        expX[i] = pow( x[i] * m_posScale, inAlpha );
        //expX[i] = expTable( x[i] / inAlpha );
        /*if( expX[i] == 0 )
    {
        printf( "ERR x[i] %f alpha %f \n", x[i], inAlpha );
    }
    assert( expX[i] != 0 );*/
    }
}

void MyNLP::UpdateExpValueForEachPin( const int& n, const double* x, double* expPins, const double& inAlpha )
{
    for( unsigned int pinId=0; pinId<m_pDB->m_pins.size(); pinId++ )
    {
        int blockId = m_pDB->m_pins[pinId].moduleId;

        // 2006-02-20
        if( m_usePin[blockId] == false )
            continue;	// save time

        double xx = x[ 2*blockId ]   + m_pDB->m_pins[ pinId ].xOff;
        double yy = x[ 2*blockId+1 ] + m_pDB->m_pins[ pinId ].yOff;
        expPins[2*pinId]   = pow( xx * m_posScale, inAlpha );
        expPins[2*pinId+1] = pow( yy * m_posScale, inAlpha );
        assert( expPins[2*pinId] != 0 );
        assert( expPins[2*pinId+1] != 0 );
        //expPins[2*pinId]   = expTable( xx / inAlpha );
        //expPins[2*pinId+1] = expTable( yy / inAlpha );
    }
}

void MyNLP::UpdateNetsSumExp( const double* x, const double* expX )
{
    double sum_exp_xi_over_alpha;
    double sum_exp_inv_xi_over_alpha;
    double sum_exp_yi_over_alpha;
    double sum_exp_inv_yi_over_alpha;
    for( unsigned int n=0; n<m_pDB->m_nets.size(); n++ )
    {
        if( m_pDB->m_nets[n].size() == 0 )
            continue;
        calc_sum_exp_using_pin(
                    m_pDB->m_nets[n].begin(), m_pDB->m_nets[n].end(), x, expX,
                    sum_exp_xi_over_alpha, sum_exp_inv_xi_over_alpha,
                    sum_exp_yi_over_alpha, sum_exp_inv_yi_over_alpha );

        m_nets_sum_exp_xi_over_alpha[n]     = sum_exp_xi_over_alpha;
        m_nets_sum_exp_yi_over_alpha[n]     = sum_exp_yi_over_alpha;
        m_nets_sum_exp_inv_xi_over_alpha[n] = sum_exp_inv_xi_over_alpha;
        m_nets_sum_exp_inv_yi_over_alpha[n] = sum_exp_inv_yi_over_alpha;
    }

    for( unsigned int n=0; n<m_pDB->m_nets.size(); n++ )
    {
        if( m_pDB->m_nets[n].size() == 0 )
            continue;
        m_nets_sum_p_x_pos[n]     = pow( m_nets_sum_exp_xi_over_alpha[n], 1/_alpha-1 );
        m_nets_sum_p_inv_x_neg[n] = pow( m_nets_sum_exp_inv_xi_over_alpha[n], -1/_alpha-1 );
        m_nets_sum_p_y_pos[n]     = pow( m_nets_sum_exp_yi_over_alpha[n], 1/_alpha-1 );
        m_nets_sum_p_inv_y_neg[n] = pow( m_nets_sum_exp_inv_yi_over_alpha[n], -1/_alpha-1 );

        // not used --hdgao
        m_nets_sum_p_inv_x_pos[n] = pow( m_nets_sum_exp_inv_xi_over_alpha[n], 1/_alpha-1 );
        m_nets_sum_p_inv_y_pos[n] = pow( m_nets_sum_exp_inv_yi_over_alpha[n], 1/_alpha-1 );
        m_nets_sum_p_x_neg[n]     = pow( m_nets_sum_exp_xi_over_alpha[n], -1/_alpha-1 );
        m_nets_sum_p_y_neg[n]     = pow( m_nets_sum_exp_yi_over_alpha[n], -1/_alpha-1 );


    }

}

double MyNLP::GetWL( const int& n, const double* x, const double* expX, const double& alpha )
{
    totalWL = 0;
    for( unsigned int n=0; n<m_pDB->m_nets.size(); n++ )	// for each net
    {
        if( m_pDB->m_nets[n].size() == 0 )
            continue;

        double invAlpha = 1.0 / alpha;
        totalWL +=
                pow( m_nets_sum_exp_xi_over_alpha[n], invAlpha ) -
                pow( m_nets_sum_exp_inv_xi_over_alpha[n], -invAlpha ) +
                pow( m_nets_sum_exp_yi_over_alpha[n], invAlpha ) -
                pow( m_nets_sum_exp_inv_yi_over_alpha[n], -invAlpha );
    }
    return totalWL / m_posScale;
}


bool MyNLP::eval_f(int n, const double* x, const double* expX, bool new_x, double& obj_value)
{
    double time_start = seconds();
    
    totalWL = GetWL( n, x, expX, _alpha );
    time_wl += seconds() - time_start;
    
    double time_start_2 = seconds();
    density = GetDensityPanelty();
    time_update_grid += seconds() - time_start_2;

    obj_value = (totalWL * _weightWire) + (density * _weightDensity);
    
    time_f += seconds() - time_start;
    return true;
}

bool MyNLP::eval_f_HPWL(int n, const double* x, const double* expX, bool new_x, double& obj_value)
{
    double time_start = seconds();
    
    UpdateBlockPosition( x );
    totalWL = m_pDB->CalcHPWL();
    time_wl += seconds() - time_start;
    
    double time_start_2 = seconds();
    density = GetDensityPanelty();
    time_update_grid += seconds() - time_start_2;

    obj_value = (totalWL * _weightWire) + (density * _weightDensity);
    
    time_f += seconds() - time_start;
    return true;
}

void MyNLP::PrintPotentialGrid()
{
    for( int i=(int)m_gridPotential.size()-1; i>=0; i-- )
    {
        for( unsigned int j=0; j<m_gridPotential[i].size(); j++ )
        {
            printf( "%4.1f ", (m_gridPotential[i][j]-m_expBinPotential[i][j])/m_expBinPotential[i][j] );
        }
        printf( "\n" );
    }
    printf( "\n\n" );
}


double MyNLP::GetDensityPanelty()
{
    double density = 0;
    for( unsigned int i=0; i<m_gridPotential.size(); i++ )
    {
        for( unsigned int j=0; j<m_gridPotential[i].size(); j++ )
        {
            density += ( m_gridPotential[i][j] - m_expBinPotential[i][j] ) *
                    ( m_gridPotential[i][j] - m_expBinPotential[i][j] );
        }
    }
    return density;
}

void MyNLP::InitModuleNetPinId()
{
    //printf( "Init module-net-pin id\n" );
    m_moduleNetPinId.resize( m_pDB->m_modules.size() );
    for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
    {
        m_moduleNetPinId[i].resize( m_pDB->m_modules[i].m_netsId.size() );
        for( unsigned int j=0; j<m_pDB->m_modules[i].m_netsId.size(); j++ )
        {
            int netId = m_pDB->m_modules[i].m_netsId[j];
            int pinId = -1;
            for( unsigned int p=0; p<m_pDB->m_nets[netId].size(); p++ )
            {
                if( m_pDB->m_pins[ m_pDB->m_nets[netId][p] ].moduleId == (int)i )
                {
                    pinId = m_pDB->m_nets[netId][p];
                    break;
                }
            }
            assert( pinId != -1 );  // floating pin? (impossible for bookshelf format)
            m_moduleNetPinId[i][j] = pinId;
        } // each net to the module
    } // each module
}

bool MyNLP::eval_grad_f(int n, const double* x, const double* expX, bool new_x, double* grad_f)
{
    double time_used = seconds();

    // grad WL
    if( _weightWire > 0 )	//TEST
        for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )	// for each block
        {
            if( m_pDB->m_modules[i].m_isFixed || m_pDB->m_modules[i].m_netsId.size() == 0 )
                continue;


            grad_wire[ 2*i ] = 0;
            grad_wire[ 2*i+1 ] = 0;

            for( unsigned int j=0; j<m_pDB->m_modules[i].m_netsId.size(); j++ )
            {
                // for each net connecting to the block
                size_t netId = m_pDB->m_modules[i].m_netsId[j];
                if( m_pDB->m_nets[netId].size() == 0 ) // floating-module
                    continue;

                // TODO: modification for LEF/DEF input
                // no floating pin for bookshelf format
                //if( m_pDB->m_nets[netId].size() == 0 )
                //	continue;

                int selfPinId = m_moduleNetPinId[i][j];

                if( m_usePin[i] )
                {
                    assert( selfPinId != -1 );
                    double xx = x[ 2*i ]   + m_pDB->m_pins[ selfPinId ].xOff;
                    double yy = x[ 2*i+1 ] + m_pDB->m_pins[ selfPinId ].yOff;
                    xx *= m_posScale;
                    yy *= m_posScale;

                    grad_wire[ 2*i ] +=
                            m_nets_sum_p_x_pos[netId]     * _expPins[2*selfPinId] / xx -
                            m_nets_sum_p_inv_x_neg[netId] / _expPins[2*selfPinId] / xx;
                    grad_wire[ 2*i+1 ] +=
                            m_nets_sum_p_y_pos[netId]     * _expPins[2*selfPinId+1] / yy -
                            m_nets_sum_p_inv_y_neg[netId] / _expPins[2*selfPinId+1] / yy;
                }
                else
                {
                    double xx = x[ 2*i ];
                    double yy = x[ 2*i+1 ];
                    xx *= m_posScale;
                    yy *= m_posScale;

                    grad_wire[ 2*i ] +=
                            m_nets_sum_p_x_pos[netId]     * _expX[2*i] / xx  -
                            m_nets_sum_p_inv_x_neg[netId] / _expX[2*i] / xx;
                    grad_wire[ 2*i+1 ] +=
                            m_nets_sum_p_y_pos[netId]     * _expX[2*i+1] / yy -
                            m_nets_sum_p_inv_y_neg[netId] / _expX[2*i+1] / yy;
                }

            } // for each pin in the module
        } // for each module
    time_grad_wl += seconds() - time_used;



    // grad Density
    double time_start_2 = seconds();
#if 0
    UpdateBinForce();	// diffusion
#endif
    
    double gradDensityX;
    double gradDensityY;
    for( int i=0; i<(int)m_pDB->m_modules.size(); i++ )	    // for each cell
    {
        if( m_pDB->m_modules[i].m_isFixed )
            continue;
#if 0
        GetDiffusionGrad( x, i, gradDensityX, gradDensityY );	    // diffusion
        grad_f[2*i]   = -(2 * gradDensityX);
        grad_f[2*i+1] = -(2 * gradDensityY);
#endif
#if 1	
        GetPotentialGrad( x, i, gradDensityX, gradDensityY );	    // bell-shaped potential
        grad_potential[2*i]   = /*2 **/ gradDensityX;
        grad_potential[2*i+1] = /*2 **/ gradDensityY;
#endif	
    } // for each cell
    time_grad_potential += seconds() - time_start_2;
    
    // compute total fouce
    for( int i =0; i<n; i++ )
        grad_f[i] = _weightDensity * grad_potential[i] + grad_wire[i] * _weightWire;
    
    time_grad_f += seconds()-time_used;
    return true;
}

/*
void MyNLP::GetDensityGrad( const double* x, const int& b, double& gradX, double& gradY )
{
    double w  = m_pDB->m_modules[b].m_width;
    double h  = m_pDB->m_modules[b].m_height;

    // bottom-left
    double left   = x[b*2]   - w * 0.5;
    double bottom = x[b*2+1] - h * 0.5;
    double right  = left   + w;
    double top    = bottom + h;

    // find nearest gird
    int gLeft   = static_cast<int>( floor( (left - m_pDB->m_coreRgn.left) / m_gridDensityWidth ) );
    int gBottom = static_cast<int>( floor( (bottom - m_pDB->m_coreRgn.bottom) / m_gridDensityHeight ) );
    int gRight  = static_cast<int>( floor( (right - m_pDB->m_coreRgn.left) / m_gridDensityWidth ) );
    int gTop    = static_cast<int>( floor( (top - m_pDB->m_coreRgn.bottom) / m_gridDensityHeight ) );

    gradX = 0;
    gradY = 0;
    for( int xOff = gLeft; xOff < gRight; xOff++ )
    {
        for( int yOff = gBottom; yOff < gTop; yOff++ )
        {
        gradX += ( m_gridDensity[xOff+1][yOff] - m_gridDensity[xOff][yOff] ) / m_gridDensityWidth;
        gradY += ( m_gridDensity[xOff][yOff+1] - m_gridDensity[xOff][yOff] ) / m_gridDensityHeight;
        }
    }
    gradX /= m_gridDensityTarget;
    gradY /= m_gridDensityTarget;
}*/


void MyNLP::GetDiffusionGrad( const double* x, const int& i, double& gradX, double& gradY )
{
    double cellX = x[i*2];
    double cellY = x[i*2+1];

    int gx, gy;	// left bottom grid
    GetClosestGrid( cellX, cellY, gx, gy );

    double xx, yy;	// left bottom grid center coordiante (x, y)
    xx = GetXGrid( gx );
    yy = GetYGrid( gy );
    if( xx > cellX )
    {
        assert( gx > 0 );
        gx--;
        xx -= m_potentialGridWidth;
    }
    if( yy > cellY )
    {
        assert( gy > 0 );	// TODO boundary
        gy--;
        yy -= m_potentialGridHeight;
    }

    // interpolation
    double alpha = ( cellX - xx ) / m_potentialGridWidth;	// x-direction
    double beta  = ( cellY - yy ) / m_potentialGridHeight;	// y-direction
    gradX = m_binForceX[gx][gy] +
            alpha * (m_binForceX[gx+1][gy] - m_binForceX[gx][gy]) +
            beta  * (m_binForceX[gx][gy+1] - m_binForceX[gx][gy]) +
            alpha * beta * (m_binForceX[gx][gy] + m_binForceX[gx+1][gy+1] -
            m_binForceX[gx+1][gy] - m_binForceX[gx][gy+1] );
    gradY = m_binForceY[gx][gy] +
            alpha * (m_binForceY[gx+1][gy] - m_binForceY[gx][gy]) +
            beta  * (m_binForceY[gx][gy+1] - m_binForceY[gx][gy]) +
            alpha * beta * (m_binForceY[gx][gy] + m_binForceY[gx+1][gy+1] -
            m_binForceY[gx+1][gy] - m_binForceY[gx][gy+1] );
}


void MyNLP::GetPotentialGrad( const double* x, const int& i, double& gradX, double& gradY )
{
    double cellX = x[i*2];
    double cellY = x[i*2+1];

    double width  = m_pDB->m_modules[i].m_width;
    double height = m_pDB->m_modules[i].m_height;
    double left   = cellX - width  * 0.5 - _potentialRX;
    double bottom = cellY - height * 0.5 - _potentialRY;
    double right  = cellX + ( cellX - left );
    double top    = cellY + ( cellY - bottom );
    if( left   < m_pDB->m_coreRgn.left )	left   = m_pDB->m_coreRgn.left;
    if( bottom < m_pDB->m_coreRgn.bottom )	bottom = m_pDB->m_coreRgn.bottom;
    if( right  > m_pDB->m_coreRgn.right )	right  = m_pDB->m_coreRgn.right;
    if( top    > m_pDB->m_coreRgn.top )	    top    = m_pDB->m_coreRgn.top;
    int gx, gy;
    GetClosestGrid( left, bottom, gx, gy );
    assert(gx >= 0);
    assert(gy <= m_potentialGridSize);

    if( gx < 0 )	gx = 0;
    if( gy < 0 )	gy = 0;

    int gxx, gyy;
    double xx, yy;
    gradX = 0.0;
    gradY = 0.0;

    //// TEST (std-cell)
    if( height < m_potentialGridHeight && width < m_potentialGridWidth )
        width = height = 0;

    for( gxx = gx, xx = GetXGrid( gx );
         xx <= right && gxx < (int)m_gridPotential.size();
         gxx++, xx += m_potentialGridWidth )
    {

        for( gyy = gy, yy = GetYGrid( gy );
             yy <= top && gyy < (int)m_gridPotential.size() ;
             gyy++, yy += m_potentialGridHeight )
        {

            double gX = 0;
            double gY = 0;
            // TEST
            //if( m_gridPotential[ gxx ][ gyy ] > m_expBinPotential[gxx][gyy] )  // TEST for ispd05
            {
                gX = ( m_gridPotential[gxx][gyy] - m_expBinPotential[gxx][gyy] ) *
                        _cellPotentialNorm[i] *
                        GetGradPotential( cellX, xx, _potentialRX, width ) *
                        GetPotential(     cellY, yy, _potentialRY, height );
                gY =  ( m_gridPotential[gxx][gyy] - m_expBinPotential[gxx][gyy] ) *
                        _cellPotentialNorm[i] *
                        GetPotential(     cellX, xx, _potentialRX, width  ) *
                        GetGradPotential( cellY, yy, _potentialRY, height );
            }

            gradX += gX;
            gradY += gY;
        }
    } // for each grid
}


void MyNLP::calc_sum_exp_using_pin( 
        const vector<int>::const_iterator& begin, const vector<int>::const_iterator& end,
        const double* x, const double* expX,
        double& sum_exp_xi_over_alpha, double& sum_exp_inv_xi_over_alpha,
        double& sum_exp_yi_over_alpha, double& sum_exp_inv_yi_over_alpha, int id )
{
    double t_start = seconds();
    
    sum_exp_xi_over_alpha = 0;
    sum_exp_inv_xi_over_alpha = 0;
    sum_exp_yi_over_alpha = 0;
    sum_exp_inv_yi_over_alpha = 0;

    vector<int>::const_iterator ite;
    int pinId;
    int blockId;
    for( ite=begin; ite!=end; ++ite )
    {
        // for each pin of the net
        pinId   = *ite;
        blockId = m_pDB->m_pins[ pinId ].moduleId;

        /*sum_exp_xi_over_alpha     += expX[2*blockId];
    sum_exp_inv_xi_over_alpha += 1.0 / expX[2*blockId];
    sum_exp_yi_over_alpha     += expX[2*blockId+1];
    sum_exp_inv_yi_over_alpha += 1.0 / expX[2*blockId+1];
    */
#if 1
        if( m_usePin[blockId] /*&& blockId != id*/ )	// macro or self pin
            //if( blockId != id )
        {
            // handle pins
            sum_exp_xi_over_alpha     += _expPins[ 2*pinId ];
            sum_exp_inv_xi_over_alpha += 1.0 / _expPins[ 2*pinId ];
            sum_exp_yi_over_alpha     += _expPins[ 2*pinId+1 ];
            sum_exp_inv_yi_over_alpha += 1.0 / _expPins[ 2*pinId+1 ];
        }
        else
        {
            // use block center
            //assert( expX[2*blockId] != 0);
            //assert( expX[2*blockId+1] != 0 );
            sum_exp_xi_over_alpha     += expX[2*blockId];
            sum_exp_inv_xi_over_alpha += 1.0 / expX[2*blockId];
            sum_exp_yi_over_alpha     += expX[2*blockId+1];
            sum_exp_inv_yi_over_alpha += 1.0 / expX[2*blockId+1];
        }
#endif
    }
    time_sum_exp += seconds() - t_start;
} 
/*
void MyNLP::finalize_solution(SolverReturn status,
                              Index n, const double* x, const double* z_L, const double* z_U,
                              Index m, const_potentialRX double* g, const double* lambda,
                              double obj_value)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution. Since the solution is displayed to the console,
  // we currently do nothing here.
  

    UpdateBlockPosition( x );
    printf( "HPWL= %f\n", m_pDB->CalcHPWL() );
    m_pDB->OutputGnuplotFigure( "fig_final.plt", false );

}*/

void MyNLP::UpdateBlockPosition( const double* x )
{
    for( int i=0; i<(int)m_pDB->m_modules.size(); i++ )
    {
        if( m_pDB->m_modules[i].m_isFixed == false )
        {
            m_pDB->MoveModuleCenter( i, x[i*2], x[i*2+1] );
        }
    }
}

void MyNLP::CreatePotentialGrid()
{
    //printf( "Create Potential Grid\n" );
    m_gridPotential.clear(); // remove old values
    
    int realGridSize = m_potentialGridSize;

    m_gridPotential.resize( realGridSize );
    m_basePotential.resize( realGridSize );
    for( unsigned int i=0; i<m_gridPotential.size(); i++ )
    {
        m_basePotential[i].resize( realGridSize, 0 );
        m_gridPotential[i].resize( realGridSize, 0 );
    }
    
    m_potentialGridWidth  = ( m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left ) / m_potentialGridSize;
    m_potentialGridHeight = ( m_pDB->m_coreRgn.top   - m_pDB->m_coreRgn.bottom ) / m_potentialGridSize;
    _potentialRX = m_potentialGridWidth  * _potentialGridR;
    _potentialRY = m_potentialGridHeight * _potentialGridR;

}


void MyNLP::ClearPotentialGrid()
{
    for( int gx=0; gx<(int)m_gridPotential.size(); gx++ )
        fill( m_gridPotential[gx].begin(), m_gridPotential[gx].end(), 0.0 );
}

#if 0
void MyNLP::UpdateBinForce()	// 2006-02-21
{
    for( unsigned int i=0; i<m_gridPotential.size(); i++ )
        for( unsigned int j=0; j<m_gridPotential[i].size(); j++ )
        {
            if( j == 0 || j == m_gridPotential.size()-1 )    // left and right
                m_binForceY[i][j] = 0;
            else
                m_binForceY[i][j] = - ( m_gridPotential[i][j+1] - m_gridPotential[i][j-1] ) /
                        m_gridPotential[i][j] / 2;

            if( i == 0 || i == m_gridPotential.size()-1 )    // bottom and top
                m_binForceX[i][j] = 0;
            else
                m_binForceX[i][j] = - ( m_gridPotential[i+1][j] - m_gridPotential[i-1][j] ) /
                        m_gridPotential[i][j] / 2;
        }
}
#endif


/*
void MyNLP::UpdateGridPotentialByGrad()
    //double xx = x[ 2*blockId ]   + m_pDB->m_pins[ pinId ].xOff;
    //double yy = x[ 2*blockId+1 ] + m_pDB->m_pins[ pinId ].yOff;
{
    for( unsigned int i=0; i<m_gridPotential.size(); i++ )
    for( unsigned int j=0; j<m_gridPotential[i].size(); j++ )
        m_gridPotential[i][j] -= _gridGradPotential[i][j] * stepSize * _weightDensity;
}
*/

void MyNLP::UpdateExpBinPotential( double util )
{
    double binArea = m_potentialGridWidth * m_potentialGridHeight;

    if( util < 0 )
        util = 1.0; // use all space

    double totalFree = 0;
    int zeroSpaceBin = 0;
    m_expBinPotential.resize( m_basePotential.size() );
    for( unsigned int i=0; i<m_basePotential.size(); i++ )
    {
        m_expBinPotential[i].resize( m_basePotential[i].size() );
        for( unsigned int j=0; j<m_basePotential[i].size(); j++ )
        {
            double base = m_basePotential[i][j];
            double free = binArea - base;
            if( free > 1e-4 )
            {
                m_expBinPotential[i][j] = free * util;
                totalFree += m_expBinPotential[i][j];
            }
            else
            {
                m_expBinPotential[i][j] = 0.0;
                zeroSpaceBin++;
            }
        }
    }

    if( param.bShow )
    {
        printf( "PBIN: Expect bin potential utilization: %f\n", util );
        printf( "PBIN: Zero space bin # = %d\n", zeroSpaceBin );
        printf( "PBIN: Total free potential = %.0f (%.5f)\n", totalFree, m_pDB->m_totalMovableModuleArea / totalFree );
    }

    // TODO: scaling?
    //assert( m_pDB->m_totalMovableModuleArea / totalFree <= 1.000001 );
    double alwaysOver = 0.0;
    if( m_targetUtil > 0.0 && m_targetUtil < 1.0 )
    {
        for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
        {
            if( m_pDB->m_modules[i].m_isFixed )
                continue;
            if( m_pDB->m_modules[i].m_width >= 2 * m_potentialGridWidth &&
                    m_pDB->m_modules[i].m_height >= 2 * m_potentialGridHeight )
            {
                alwaysOver +=
                        (m_pDB->m_modules[i].m_width - m_potentialGridWidth ) *
                        (m_pDB->m_modules[i].m_height - m_potentialGridHeight ) *
                        (1.0 - m_targetUtil );
            }
        }
        if( param.bShow )
            printf( "PBIN: Always over: %.0f (%.1f%%)\n", alwaysOver, alwaysOver/m_pDB->m_totalMovableModuleArea*100.0 );
    }
    m_alwaysOverPotential = alwaysOver;
}

void MyNLP::SmoothPotentialBase( const double& delta )
{
    
    // find the max potential (TODO: comnpute one time is enough)
    double maxPotential = 0;
    double avgPotential = 0;
    double totalPotential = 0;
    for( unsigned int i=0; i<m_basePotentialOri.size(); i++ )
        for( unsigned int j=0; j<m_basePotentialOri[i].size(); j++ )
        {
            totalPotential += m_basePotentialOri[i][j];
            if( m_basePotentialOri[i][j] > maxPotential )
                maxPotential = m_basePotentialOri[i][j];
        }
    avgPotential = totalPotential / (m_basePotentialOri.size() * m_basePotentialOri.size() );

    if( totalPotential == 0 )
        return; // no preplaced
    
    // apply TSP-style smoothing
    double newTotalPotential = 0;
    for( unsigned int i=0; i<m_basePotential.size(); i++ )
        for( unsigned int j=0; j<m_basePotential[i].size(); j++ )
        {
            if( m_basePotentialOri[i][j] >= avgPotential )
            {
                m_basePotential[i][j] =
                        avgPotential +
                        pow( ( m_basePotentialOri[i][j] - avgPotential ) / maxPotential, delta ) * maxPotential;
            }
            else
            {
                m_basePotential[i][j] =
                        avgPotential -
                        pow( ( avgPotential - m_basePotentialOri[i][j] ) / maxPotential, delta ) * maxPotential;
            }
            newTotalPotential += m_basePotential[i][j];
        }
    
    // normalization
    double ratio = totalPotential / newTotalPotential;
    for( unsigned int i=0; i<m_basePotential.size(); i++ )
        for( unsigned int j=0; j<m_basePotential[i].size(); j++ )
            m_basePotential[i][j] = m_basePotential[i][j] * ratio;

    //printf( "Smooth %.0f (%.0f->%.0f)\n", delta, totalPotential, newTotalPotential );
}

void MyNLP::UpdatePotentialGridBase( const double* x )
{
    double time_start = seconds();

    double binArea = m_potentialGridWidth * m_potentialGridHeight;
    m_binFreeSpace.resize( m_basePotential.size() );
    for( unsigned int i=0; i<m_basePotential.size(); i++ )
    {
        fill( m_basePotential[i].begin(), m_basePotential[i].end(), 0.0 );
        m_binFreeSpace[i].resize( m_basePotential[i].size() );
        fill( m_binFreeSpace[i].begin(), m_binFreeSpace[i].end(), binArea );
    }

    for( int i=0; i<(int)m_pDB->m_modules.size(); i++ )
    {
        // for each cell. cell ci coordinate is ( x[i*2], x[i*2+1] )

        if( m_pDB->m_modules[i].m_isFixed == false )
            continue;

        // TODO: BUG when shrinking core?
        if( m_pDB->m_modules[i].m_isOutCore )
            continue;	// pads?

        int gx, gy;
        double cellX = x[i*2];
        double cellY = x[i*2+1];
        double width  = m_pDB->m_modules[i].m_width;
        double height = m_pDB->m_modules[i].m_height;

        double potentialRX = _potentialRX;
        double potentialRY = _potentialRY;
        //double left   = cellX - width * 0.5  - potentialRX;
        //double bottom = cellY - height * 0.5 - potentialRY;
        double left   = cellX - width * 0.5;  // for gaussian smoothing
        double bottom = cellY - height * 0.5; // for gaussian smoothing
        double right  = cellX + (cellX - left);
        double top    = cellY + (cellY - bottom);
        if( left   < m_pDB->m_coreRgn.left )     left   = m_pDB->m_coreRgn.left;
        if( bottom < m_pDB->m_coreRgn.bottom )   bottom = m_pDB->m_coreRgn.bottom;
        if( top    > m_pDB->m_coreRgn.top )      top    = m_pDB->m_coreRgn.top;
        if( right  > m_pDB->m_coreRgn.right )    right  = m_pDB->m_coreRgn.right;
        GetClosestGrid( left, bottom, gx, gy );
        if( gx < 0 )  gx = 0;
        if( gy < 0 )  gy = 0;

        double totalPotential = 0;
        vector< potentialStruct > potentialList;
        int gxx, gyy;
        double xx, yy;

        //if( m_useBellPotentialForPreplaced == false )
        {
            // "Exact density for the potential"
            for( gxx = gx, xx = GetXGrid(gx); xx<=right ; gxx++, xx+=m_potentialGridWidth )
            {
                for( gyy = gy, yy = GetYGrid(gy); yy<=top ; gyy++, yy+=m_potentialGridHeight )
                {
                    m_basePotential[gxx][gyy] +=
                            getOverlap( left, right, xx, xx+m_potentialGridWidth ) *
                            getOverlap( bottom, top, yy, yy+m_potentialGridHeight );

                    m_binFreeSpace[gxx][gyy] -=
                            getOverlap( left, right, xx, xx+m_potentialGridWidth ) *
                            getOverlap( bottom, top, yy, yy+m_potentialGridHeight );
                }
            }
            continue;
        }

        for( gxx = gx, xx = GetXGrid(gx); xx<=right ; gxx++, xx+=m_potentialGridWidth )
        {
            for( gyy = gy, yy = GetYGrid(gy); yy<=top ; gyy++, yy+=m_potentialGridHeight )
            {
                double potential = GetPotential( cellX, xx, potentialRX, width ) *
                        GetPotential( cellY, yy, potentialRY, height );
                if( potential > 0 )
                {
                    totalPotential += potential;
                    potentialList.push_back( potentialStruct( gxx, gyy, potential ) );
                }
            }
        }

        // normalize the potential so that total potential equals the cell area
        double scale = m_pDB->m_modules[i].m_area / totalPotential;
        //printf( "totalPotential = %f\n", totalPotential );

        _cellPotentialNorm[i] = scale;	    // normalization factor for the cell i

        vector< potentialStruct >::const_iterator ite;
        for( ite=potentialList.begin(); ite!=potentialList.end(); ++ite )
        {
            if(	ite->gx < 0 || ite->gx >= (int)m_gridPotential.size() ||
                    ite->gy < 0 || ite->gy >= (int)m_gridPotential[ite->gx].size() )
                continue; // bin may be outside when core-shrinking is applied
            else
                m_basePotential[ ite->gx ][ ite->gy ] += ite->potential * scale;
        }


    } // for each cell
    time_up_potential += seconds() - time_start;

    m_basePotentialOri = m_basePotential;   // make a copy for TSP-style smoothing

    /*
    double totalFreeSpace = 0;
    for( unsigned int i=0; i<m_binFreeSpace.size(); i++ )
    {
    for( unsigned int j=0; j<m_binFreeSpace[i].size(); j++ )
    {
        if( m_binFreeSpace[i][j] < 0 )
        m_binFreeSpace[i][j] = 0;
        totalFreeSpace += m_binFreeSpace[i][j];
    }
    }
    printf( "totalFreeSpace: %.0f\n", totalFreeSpace );
    */
}


void MyNLP::UpdatePotentialGrid( const double* x )
{
    double time_start = seconds();
    ClearPotentialGrid();
    for( int i=0; i<(int)m_pDB->m_modules.size(); i++ )
    {
        // for each cell. cell ci coordinate is ( x[i*2], x[i*2+1] )

        if( m_pDB->m_modules[i].m_isOutCore )
            continue;

        // preplaced blocks are stored in m_basePotential
        if( m_pDB->m_modules[i].m_isFixed )
            continue;

        int gx, gy;
        double cellX = x[i*2];
        double cellY = x[i*2+1];
        double potentialRX = _potentialRX;// 2 * m_PotentialGridWidth;
        double potentialRY = _potentialRY;// 2 * m_PotentialGridHeight
        double width  = m_pDB->m_modules[i].m_width;
        double height = m_pDB->m_modules[i].m_height;
        double left   = cellX - width * 0.5  - potentialRX;
        double bottom = cellY - height * 0.5 - potentialRY;
        double right  = cellX + (cellX - left);
        double top    = cellY + (cellY - bottom);
        if( left   < m_pDB->m_coreRgn.left )     left   = m_pDB->m_coreRgn.left;
        if( bottom < m_pDB->m_coreRgn.bottom )   bottom = m_pDB->m_coreRgn.bottom;
        if( top    > m_pDB->m_coreRgn.top )      top    = m_pDB->m_coreRgn.top;
        if( right  > m_pDB->m_coreRgn.right )    right  = m_pDB->m_coreRgn.right;
        GetClosestGrid( left, bottom, gx, gy );

        double totalPotential = 0;
        vector< potentialStruct > potentialList;
        int gxx, gyy;
        double xx, yy;

        //// TEST (convert to std-cell)
        if( height < m_potentialGridHeight && width < m_potentialGridWidth )
            width = height = 0;
        // xx,yy the center position of bin(gx,gy)
        for( gxx = gx, xx = GetXGrid(gx); xx<=right ; gxx++, xx+=m_potentialGridWidth )
        {// m_pDB->m_coreRgn.left + gx * m_potentialGridWidth + 0.5 * m_potentialGridWidth;
            for( gyy = gy, yy = GetYGrid(gy); yy<=top ; gyy++, yy+=m_potentialGridHeight )
            {
                double potential = GetPotential( cellX, xx, potentialRX, width ) *
                        GetPotential( cellY, yy, potentialRY, height );
                if( potential > 0 )
                {
                    totalPotential += potential;
                    potentialList.push_back( potentialStruct( gxx, gyy, potential ) );
                }
#if DEBUG_hd
                printf("gx = %-4d,gy = %-4d,potential = %-4.3f\n",gxx,gyy,potential);
                fflush(stdout);
#endif
            }
        }

        // normalize the potential so that total potential equals the cell area
        double scale = m_pDB->m_modules[i].m_area / totalPotential;
        //printf( "totalPotential = %f\n", totalPotential );

        _cellPotentialNorm[i] = scale;	    // normalization factor for the cell i
        vector< potentialStruct >::const_iterator ite;
        for( ite=potentialList.begin(); ite!=potentialList.end(); ++ite )
        {
#if 0	    
            assert( ite->gx <=  (int)m_gridPotential.size() );
            assert( ite->gy <=  (int)m_gridPotential.size() );
            assert( ite->gx >= 0 );
            assert( ite->gy >= 0 );
#endif
            m_gridPotential[ ite->gx ][ ite->gy ] += ite->potential * scale;
        }

    } // for each cell
    time_up_potential += seconds() - time_start;

}

/*double MyNLP::GetGridWidth()
{
    return m_potentialGridWidth;
}*/

/*double MyNLP::GetPotential( const double& x1, const double& x2, const double& r )
{
    double d = fabs( x1 - x2 );

    if( d <= r * 0.5 )
    return 1.0 - 2 * d * d / ( r * r );
    else if( d <= r )
    return 2 * ( d - r ) * ( d - r ) / ( r * r );
    else
    return 0;
}*/

double MyNLP::GetPotential( const double& x1, const double& x2, const double& r, const double& w )
{
    double d = fabs( x1 - x2 );
    double a = 4.0 / ( w + r ) / ( w + 2 * r );
    double b = 4.0 / r / ( w + 2.0 * r );
    
    if( d <= w * 0.5 + r * 0.5 )
        return 1.0 - a * d * d;
    else if( d <= w * 0.5 + r )
        return b * ( d - r - w * 0.5 ) * ( d - r - w * 0.5);
    else
        return 0.0;
}
/*
double MyNLP::GetGradPotential( const double& x1, const double& x2, const double& r )
{
    double d;
    if( x1 >= x2 )  // right half
    {
    d = x1 - x2;	// d >= 0
    if( d <= r * 0.5 )
        return -4.0 * d / ( r * r );
    else if( d <= r )
        return +4.0 * ( d - r ) / ( r * r );
    else
        return 0;
    }
    else    // left half
    {
    d = x2 - x1;	// d >= 0
    if( d <= r * 0.5 )
        return +4.0 * d / ( r * r );
    else if( d <= r )
        return -4.0 * ( d - r ) / ( r * r );
    else
        return 0;
    }
}*/

double MyNLP::GetGradPotential( const double& x1, const double& x2, const double& r, const double& w )
{
    //double w = 0;
    double d;
    double a = 4.0 / ( w + r ) / ( w + 2.0 * r );
    double b = 4.0 / r / ( w + 2.0 * r );

    if( x1 >= x2 )  // right half
    {
        d = x1 - x2;	// d >= 0
        if( d <= w * 0.5 + r * 0.5 )
            return -2.0 * a * d;
        else if( d <= w * 0.5 + r )
            return +2.0 * b * ( d - r - w * 0.5);
        else
            return 0;
    }
    else    // left half
    {
        d = x2 - x1;	// d >= 0
        if( d <= w * 0.5 + r * 0.5 )
            return +2.0 * a * d;
        else if( d <= w * 0.5 + r )
            return -2.0 * b * ( d - r - w * 0.5);
        else
            return 0;
    }
}

/*double MyNLP::GetGradGradPotential( const double& x1, const double& x2, const double& r )
{
    double d = fabs( x1 - x2 );

    if( d <= r * 0.5 )
    return -4.0 / ( r * r );
    else if( d <= r )
    return +4.0 / ( r * r );
    else
    return 0;
}*/

/*void   MyNLP::GetGridCenter( const int& gx, const int& gy, double& x1, double& y1 )
{
    assert( gx <= m_potentialGridSize );
    assert( gy <= m_potentialGridSize );
    assert( gx >= 0 );
    assert( gy >= 0 );
    
    x1 = m_pDB->m_coreRgn.left   + gx * m_potentialGridWidth  + 0.5 * m_potentialGridWidth;
    y1 = m_pDB->m_coreRgn.bottom + gy * m_potentialGridHeight + 0.5 * m_potentialGridHeight;
}*/

double MyNLP::GetXGrid( const int& gx )
{
    return m_pDB->m_coreRgn.left + gx * m_potentialGridWidth + 0.5 * m_potentialGridWidth;
}

double MyNLP::GetYGrid( const int& gy )
{
    return  m_pDB->m_coreRgn.bottom + gy * m_potentialGridHeight + 0.5 * m_potentialGridHeight;
}

/*double MyNLP::GetPotentialToGrid( const double& x1, const int& gx, const bool& useX )
{
    if( gx < 0 || gx >= m_potentialGridSize )
    return 0.0;
    double x2, y2;
    GetGridCenter( gx, 0, x2, y2 ); // we only use "x2"
    return GetPotential( x1, x2 );
}

double MyNLP::GetGradPotentialToGrid( const double& x1, const int& gx )
{
    if( gx < 0 || gx >= m_potentialGridSize )
    return 0.0;
    double x2, y2;
    GetGridCenter( gx, 0, x2, y2 ); // we only use "x2"
    return GetGradPotential( x1, x2 );
}*/

void MyNLP::GetClosestGrid( const double& x1, const double& y1, int& gx, int& gy ) 
{
    gx = static_cast<int>( floor( ( x1 - m_pDB->m_coreRgn.left ) / m_potentialGridWidth ) );
    gy = static_cast<int>( floor( ( y1 - m_pDB->m_coreRgn.bottom ) / m_potentialGridHeight ) );

    // DEBUG
    /*if( gy >= m_gridPotential.size() || gy < 0 )
    {
    printf( "gridHeight= %f, x1= %f, y1= %f, bottom= %f, top= %f, gy= %d\n",
        m_potentialGridHeight, x1, y1, m_pDB->m_coreRgn.bottom, m_pDB->m_coreRgn.top , gy );
    }
    if( gx >= m_gridPotential.size() || gx < 0)
    {
    printf( "gridWidth = %f, y1 = %f, x1 = %f, left = %f, right = %f, gx = %d\n",
        m_potentialGridWidth, y1, x1, m_pDB->m_coreRgn.left, m_pDB->m_coreRgn.right, gx );
    }*/
    
#if 0    
    assert( gx >= 0 );
    assert( gy >= 0 );
    assert( gx < (int)m_gridPotential.size() );
    assert( gy < (int)m_gridPotential.size() );
#endif
}

void MyNLP::ClearDensityGrid()
{
    for( unsigned int i=0; i<m_gridDensity.size(); i++ )
        for( unsigned int j=0; j<m_gridDensity[i].size(); j++ )
            m_gridDensity[i][j] = 0.0;
}


void MyNLP::UpdateDensityGridSpace( const int& n, const double* x )
{
    double allSpace = m_gridDensityWidth * m_gridDensityHeight;
    for( unsigned int i=0; i<m_gridDensity.size(); i++ )
        for( unsigned int j=0; j<m_gridDensity[i].size(); j++ )
            m_gridDensitySpace[i][j] = allSpace;

    // TEST
    //return;
    
    // for each cell b, update the corresponding bin area
    for( int b=0; b<(int)m_pDB->m_modules.size(); b++ )
    {
        if( false == m_pDB->m_modules[b].m_isFixed )
            continue;

        double w  = m_pDB->m_modules[b].m_width;
        double h  = m_pDB->m_modules[b].m_height;
        double left   = x[b*2]   - w * 0.5;
        double bottom = x[b*2+1] - h * 0.5;
        double right  = left   + w;
        double top    = bottom + h;

        if( w == 0 || h == 0 )
            continue;

        // find nearest bottom-left gird
        int gx = static_cast<int>( floor( (left   - m_pDB->m_coreRgn.left)   / m_gridDensityWidth ) );
        int gy = static_cast<int>( floor( (bottom - m_pDB->m_coreRgn.bottom) / m_gridDensityHeight ) );

        if( gx < 0 )  gx = 0;
        if( gy < 0 )  gy = 0;

        for( int xOff = gx; xOff < (int)m_gridDensity.size(); xOff++ )
        {
            double binLeft  = m_pDB->m_coreRgn.left + xOff * m_gridDensityWidth;
            double binRight = binLeft + m_gridDensityWidth;
            if( binLeft >= right )
                break;

            for( int yOff = gy; yOff < (int)m_gridDensity[xOff].size(); yOff ++ )
            {
                double binBottom = m_pDB->m_coreRgn.bottom + yOff * m_gridDensityHeight;
                double binTop    = binBottom + m_gridDensityHeight;
                if( binBottom >= top )
                    break;

                m_gridDensitySpace[xOff][yOff] -=
                        getOverlap( left, right, binLeft, binRight ) *
                        getOverlap( bottom, top, binBottom, binTop );
            }
        }

    } // each module

    int zeroSpaceCount = 0;
    m_totalFreeSpace = 0;
    for( unsigned int i=0; i<m_gridDensity.size(); i++ )
        for( unsigned int j=0; j<m_gridDensity[i].size(); j++ )
        {
            if( m_gridDensitySpace[i][j] < 1e-5 )
            {
                m_gridDensitySpace[i][j] = 0.0;
                zeroSpaceCount ++;
            }
            m_totalFreeSpace += m_gridDensitySpace[i][j];
        }
    if( param.bShow )
        printf( "DBIN: zero space bins: %d\n", zeroSpaceCount );
}


void MyNLP::UpdateDensityGrid( const int& n, const double* x )
{
    ClearDensityGrid();
    
    // for each cell b, update the corresponding bin area
    for( int b=0; b<(int)m_pDB->m_modules.size(); b++ )
    {
        if(  m_pDB->m_modules[b].m_isOutCore || m_pDB->m_modules[b].m_isFixed )
            continue;

        double w  = m_pDB->m_modules[b].m_width;
        double h  = m_pDB->m_modules[b].m_height;

        // bottom-left
        double left   = x[b*2]   - w * 0.5;
        double bottom = x[b*2+1] - h * 0.5;
        double right  = left   + w;
        double top    = bottom + h;

        // find nearest gird
        int gx = static_cast<int>( floor( (left - m_pDB->m_coreRgn.left) / m_gridDensityWidth ) );
        int gy = static_cast<int>( floor( (bottom - m_pDB->m_coreRgn.bottom) / m_gridDensityHeight ) );
        if( gx < 0 ) gx = 0;
        if( gy < 0 ) gy = 0;

        // Block is always inside the core region. Do not have to check boundary.
        //double debug_area = 0;
        for( int xOff = gx; xOff < (int)m_gridDensity.size(); xOff++ )
        {
            double binLeft = m_pDB->m_coreRgn.left + m_gridDensityWidth * xOff;
            double binRight = binLeft + m_gridDensityWidth;
            if( binLeft >= right )
                break;

            for( int yOff = gy; yOff < (int)m_gridDensity[xOff].size(); yOff++ )
            {
                double binBottom = m_pDB->m_coreRgn.bottom + m_gridDensityHeight * yOff;
                double binTop    = binBottom + m_gridDensityHeight;
                if( binBottom >= top )
                    break;

                double area =
                        getOverlap( left, right, binLeft, binRight ) *
                        getOverlap( bottom, top, binBottom, binTop );

                m_gridDensity[xOff][yOff] += area;
                //debug_area += area;
            }
        }

        // TODO: check precision
        //printf( " module %d %f %f\n", b, m_pDB->m_modules[b].m_area, debug_area );

    } // each module

    /* ??? TODO: check precision
    double totalArea = 0;
    for( unsigned int i=0; i<m_gridDensity.size(); i++ )
    for( unsigned int j=0; j<m_gridDensity[i].size(); j++ )
    {
        totalArea += m_gridDensity[i][j];
    }
    printf( "%f %f\n", totalArea, m_totalMovableModuleArea );
    assert( totalArea == m_totalMovableModuleArea );*/
}

void MyNLP::CheckDensityGrid()
{
    double totalDensity = 0;
    for( int i=0; i<(int)m_gridDensity.size(); i++ )
        for( int j=0; j<(int)m_gridDensity[i].size(); j++ )
            totalDensity += m_gridDensity[i][j];

    double totalArea = 0;
    for( int i=0; i<(int)m_pDB->m_modules.size(); i++ )
    {
        if( m_pDB->m_modules[i].m_isOutCore == false )
            totalArea += m_pDB->m_modules[i].m_area;
    }

    printf( " %f %f\n", totalDensity, totalArea );
}

void MyNLP::CreateDensityGrid( int nGrid )
{
    m_gridDensity.resize( nGrid );
    for( int i=0; i<nGrid; i++ )
        m_gridDensity[i].resize( nGrid );
    
    m_gridDensitySpace.resize( nGrid );
    for( int i=0; i<nGrid; i++ )
        m_gridDensitySpace[i].resize( nGrid );
    
    m_gridDensityWidth  = ( (double)m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left ) / nGrid;
    m_gridDensityHeight = ( (double)m_pDB->m_coreRgn.top   - m_pDB->m_coreRgn.bottom ) / nGrid;
    m_gridDensityTarget = m_pDB->m_totalModuleArea / ( nGrid * nGrid );
    
    //printf( "Density Target Area = %f\n", m_gridDensityTarget );
    //printf( "Design Density = %f\n", m_gridDensityTarget/m_gridDensityWidth/m_gridDensityHeight );
    // 2006-03-21 compute always overflow area
    
    double alwaysOver = 0.0;
    if( m_targetUtil > 0.0 && m_targetUtil < 1.0 )
    {
        for( unsigned int i=0; i<m_pDB->m_modules.size(); i++ )
        {
            if( m_pDB->m_modules[i].m_isFixed )
                continue;
            if( m_pDB->m_modules[i].m_width >= 2*m_gridDensityWidth && m_pDB->m_modules[i].m_height >= 2*m_gridDensityHeight )
                alwaysOver +=
                        (m_pDB->m_modules[i].m_width - m_gridDensityWidth ) *
                        (m_pDB->m_modules[i].m_height - m_gridDensityHeight ) *
                        (1.0 - m_targetUtil );
        }
        if( param.bShow )
            printf( "DBIN: Always over: %.0f (%.1f%%)\n", alwaysOver, alwaysOver/m_pDB->m_totalMovableModuleArea*100.0 );
    }
    m_alwaysOverArea = alwaysOver;
}

/*
double MyNLP::GetDensityGridPanelty()
{
    double den = 0;
    double p;
    for( int i=0; i<(int)m_gridDensity.size(); i++ )
    {
    for( int j=0; j<(int)m_gridDensity[i].size(); j++ )
    {
        p = ( m_gridDensity[i][j] - m_gridDensityTarget ) / m_gridDensityTarget;
        p = p*p;
        den += p;
    }
    }
    return den;
}*/


double MyNLP::GetMaxDensity()
{
    double maxUtilization = 0;
    double binArea = m_gridDensityWidth * m_gridDensityHeight;
    for( int i=0; i<(int)m_gridDensity.size(); i++ )
        for( int j=0; j<(int)m_gridDensity[i].size(); j++ )
        {
            if( m_gridDensitySpace[i][j] > 1e-5 )
            {
                //double utilization = m_gridDensity[i][j] / m_gridDensitySpace[i][j];

                double preplacedArea = binArea - m_gridDensitySpace[i][j];
                double utilization = ( m_gridDensity[i][j] + preplacedArea ) / binArea;

                // TEST
                //double utilization = m_gridDensity[i][j] / m_gridDensityWidth / m_gridDensityHeight;

                if( utilization > maxUtilization )
                    maxUtilization = utilization;
            }
        }
    return maxUtilization;
}

/*
double MyNLP::GetAvgOverDensity()
{
    //const double targetDensity = 1.0;
    double avgDensity = 0;
    int overflowCount = 0;
    
    for( unsigned int i=0; i<m_gridDensity.size(); i++ )
    for( unsigned int j=0; j<m_gridDensity.size(); j++ )
        if( m_gridDensity[i][j] > m_gridDensityTarget )
        {
        overflowCount++;
            avgDensity += m_gridDensity[i][j];
        }
    return avgDensity / overflowCount / m_gridDensityTarget;
}
*/

double MyNLP::GetTotalOverDensityLB()
{
    double over = 0;
    for( unsigned int i=0; i<m_gridDensity.size(); i++ )
        for( unsigned int j=0; j<m_gridDensity.size(); j++ )
        {
            double targetSpace = m_gridDensitySpace[i][j] * m_targetUtil;
            if( targetSpace > 1e-5 && m_gridDensity[i][j]  > targetSpace  )
                over += m_gridDensity[i][j] - targetSpace;
        }

    // TODO: remove "1.0"
    return (over -m_alwaysOverArea) / (m_pDB->m_totalMovableModuleArea) + 1.0;
}


double MyNLP::GetTotalOverDensity()
{
    double over = 0;
    for( unsigned int i=0; i<m_gridDensity.size(); i++ )
        for( unsigned int j=0; j<m_gridDensity.size(); j++ )
        {
            double targetSpace = m_gridDensitySpace[i][j] * m_targetUtil;
            if( m_gridDensity[i][j]  > targetSpace  )
                over += m_gridDensity[i][j] - targetSpace;
        }

    // TODO: remove "1.0"
    return ( over - m_alwaysOverArea) / (m_pDB->m_totalMovableModuleArea) + 1.0;
}


double MyNLP::GetTotalOverPotential()
{
    double over = 0;
    for( unsigned int i=0; i<m_gridPotential.size(); i++ )
        for( unsigned int j=0; j<m_gridPotential[i].size(); j++ )
        {
            if( m_gridPotential[i][j]  > m_expBinPotential[i][j]  )
                over += m_gridPotential[i][j] - m_expBinPotential[i][j];
        }

    // TODO: remove "1.0"
    return (over - m_alwaysOverPotential) / (m_pDB->m_totalMovableModuleArea) + 1.0;
}


double MyNLP::GetNonZeroDensityGridPercent()
{
    double nonZero = 0;
    for( int i=0; i<(int)m_gridDensity.size(); i++ )
        for( int j=0; j<(int)m_gridDensity.size(); j++ )
        {
            if( m_gridDensity[i][j] > 0 ||
                    m_gridDensitySpace[i][j] == 0
                    //|| m_gridDensitySpace[i][j] < m_potentialGridWidth * m_potentialGridHeight
                    )
                nonZero += 1.0;
        }
    return nonZero / m_gridDensity.size() / m_gridDensity.size();
}


double MyNLP::GetNonZeroGridPercent()
{
    double nonZero = 0;
    for( int i=0; i<(int)m_gridPotential.size(); i++ )
        for( int j=0; j<(int)m_gridPotential.size(); j++ )
            if( m_gridPotential[i][j] > 0 )
                nonZero += 1.0;
    return nonZero / m_gridPotential.size() / m_gridPotential.size();
}


double MyNLP::GetMaxPotential()
{
    double maxDensity = 0;

    for( unsigned int i=0; i<m_gridPotential.size(); i++ )
        for( unsigned int j=0; j<m_gridPotential.size(); j++ )
            if( m_gridPotential[i][j] > maxDensity )
                maxDensity = m_gridPotential[i][j];
    return maxDensity;
}


double MyNLP::GetAvgPotential()
{
    const double targetDensity = 1.0;
    double avgDensity = 0;
    int overflowCount = 0;
    
    for( unsigned int i=0; i<m_gridPotential.size(); i++ )
        for( unsigned int j=0; j<m_gridPotential.size(); j++ )
            if( m_gridPotential[i][j] > targetDensity )
            {
                overflowCount++;
                avgDensity += m_gridPotential[i][j];
            }
    return avgDensity / overflowCount;
}

/*
// 2006-03-01
double MyNLP::GetTotalOverPotential()
{
    double totalOver = 0;
    for( unsigned int i=0; i<m_gridPotential.size(); i++ )
    for( unsigned int j=0; j<m_gridPotential.size(); j++ )
    {
        // TODO: use different expPotential in different bins
        //double targetPotential = _expPotential[i][j] * m_targetUtil;
        double targetPotential = _expPotential * m_targetUtil;
        if( m_gridPotential[i][j] > targetPotential )
        {
            totalOver += m_gridPotential[i][j] - targetPotential;
        }
    }
    // TODO: remove 1.0
    return totalOver / m_totalMovableModuleArea + 1.0;
}
*/

// Output potential data for gnuplot
void MyNLP::OutputPotentialGrid( string filename )
{
    int stepSize = (int)m_gridPotential.size() / 100;
    if( stepSize == 0 )
        stepSize = 1;
    FILE* out = fopen( filename.c_str(), "w" );
    double binArea = m_potentialGridWidth * m_potentialGridHeight;
    for( unsigned int j=0; j<m_gridPotential.size(); j+=stepSize )
    {
        for( unsigned int i=0; i<m_gridPotential.size(); i+=stepSize )
            fprintf( out, "%.03f ", (m_gridPotential[i][j] + m_basePotential[i][j]) / binArea );
        fprintf( out, "\n" );
    }
    fprintf( out, "\n" );
    fclose( out );
}


// Output potential data for gnuplot
void MyNLP::OutputDensityGrid( string filename )
{
    int stepSize = 1;
    FILE* out = fopen( filename.c_str(), "w" );
    for( unsigned int j=0; j<m_gridDensity.size(); j+=stepSize )
    {
        for( unsigned int i=0; i<m_gridDensity.size(); i+=stepSize )
        {
            double targetSpace = m_gridDensitySpace[i][j] * m_targetUtil;
            if( m_gridDensity[i][j] > targetSpace )
            {
                // % overflow
                fprintf( out, "%.03f ", (m_gridDensity[i][j]-targetSpace) / m_pDB->m_totalMovableModuleArea * 100 );
            }
            else
            {
                fprintf( out, "%.03f ", 0.0 );
            }
        }
        fprintf( out, "\n" );
    }
    fprintf( out, "\n" );
    fclose( out );
}

/*******************************************flpeng added**************************************************************/
double MyNLP::getRandNum(double x1, double x2)
{
    assert (x2 >= x1);
    int r, maxNum;
    maxNum = (int) ((x2 -x1) * 10000);
    r = rand() % maxNum;

    return (x1 + r / 10000.0);
}

double MyNLP::sg_CalcHPWL(const double* x)
{
    double HPWL = 0;
    double maxX = -10000;
    double maxY = -10000;
    double minX = DBL_MAX;
    double minY = DBL_MAX;
    double cx, cy;
    int pID, mID;
    for (unsigned int i = 0; i < m_pDB->m_nets.size(); i++)
    {
        if (m_pDB->m_nets[i].size() == 0)
            continue;

        for (unsigned int j = 0; j < m_pDB->m_nets[i].size(); j++)
        {
            pID = m_pDB->m_nets[i][j];
            mID = m_pDB->m_pins[pID].moduleId;
            assert( pID < m_pDB->m_pins.size() );
            assert( mID < m_pDB->m_modules.size() );
            cx = x[ 2 * mID ];
            cy = x[ 2 * mID + 1 ];

            minX = min( minX, cx);
            maxX = max( maxX, cx);
            minY = min( minY, cy);
            maxY = max( maxY, cy);
        }

        assert(maxX >= minX);
        assert(maxY >= minY);

        HPWL += ( (maxX - minX) + (maxY - minY) );
    }

    return HPWL;

}

void MyNLP::sg_findBeta(const int& n, const double* grad_f, const double* last_grad_f,const double* real_last_grad_f, double& beta )
{
    double l2norm = 0;
    for( int i=0; i<n; i++ )
        l2norm += real_last_grad_f[i] * real_last_grad_f[i];

    double product = 0;
    for( int i=0; i<n; i++ )
        product += grad_f[i] * ( grad_f[i] - real_last_grad_f[i] );	// g_k^T ( g_k - g_{k-1} )
    beta = product / l2norm;
}

bool MyNLP::sg_eval_f(const double* x, double& objValue)
{
    double time_start =seconds();

    //    UpdateBlockPosition(x);
    //    totalWL = m_pDB->CalcHPWL();  //CalcHPWL(x);
    double totalWL = sg_CalcHPWL(x);
    totalWL /= 1.0;

    time_wl = seconds() - time_start;

    double time_start2 = seconds();
    density = GetDensityPanelty();

    time_update_grid += seconds() - time_start2;

    objValue = totalWL * _weightWire + density * _weightDensity;

    time_f += seconds() - time_start;

    return true;
}

void MyNLP::sg_eval_grad_f(int n, const double* x, double* sub_grad_f)
{
    double time_used = seconds();
    for ( unsigned int i = 0; i < m_pDB->m_modules.size(); i++ )
    {
        sub_grad_wire[2 * i] = 0;
        sub_grad_wire[2 * i + 1] = 0;
    }

    if( _weightWire > 0 )	//TEST
        for (unsigned int i = 0; i < m_pDB->m_nets.size(); i++)
        {
            if (m_pDB->m_nets[i].size() == 0) continue;
            if (m_pDB->m_nets[i].size() == 1)
            {
                //  sg_updateSubGradWire_netBig_s(i, x);
                continue;
            }

            if (m_pDB->m_nets[i].size() == 2)
            {
                sg_updateSubGradWire_net2_s(i, x);
                // cnt2++;
                continue;
            }

            if (m_pDB->m_nets[i].size() == 3)
            {
                sg_updateSubGradWire_netBig_s(i, x);
                // cnt3++;
                continue;
            }

            if (m_pDB->m_nets[i].size() > 3)
            {
                sg_updateSubGradWire_netBig_s(i, x);
                //cnt4++;
                continue;
            }
        }

    time_grad_wl += seconds() - time_used;

    // grad Density
    double time_start_2 = seconds();

    double subGradDensityX;
    double subGradDensityY;
    for( int i=0; i<(int)m_pDB->m_modules.size(); i++ )	    // for each cell
    {
        if( m_pDB->m_modules[i].m_isFixed == true)
            continue;

        // sg_getSubGradPotential( x, i, subGradDensityX, subGradDensityY );	    // bell-shaped potential
        GetPotentialGrad( x, i, subGradDensityX, subGradDensityY );
        sub_grad_potential[2 * i]     = subGradDensityX;
        sub_grad_potential[2 * i + 1] = subGradDensityY;
    } // for each cell
    time_grad_potential += seconds() - time_start_2;

    // compute total fouce
    for( int i = 0; i < n; i++ )
        sub_grad_f[i] = _weightDensity * sub_grad_potential[i] + sub_grad_wire[i] * _weightWire;

    time_grad_f += seconds()-time_used;

}

void MyNLP::sg_updateSubGradWire_net2_s(const int &netID, const double *x,bool isPart)
{
    int pinID1 = m_pDB->m_nets[netID][0];
    int pinID2 = m_pDB->m_nets[netID][1];
    int moduleID1 = m_pDB->m_pins[pinID1].moduleId;
    int moduleID2 = m_pDB->m_pins[pinID2].moduleId;

    double xx1, xx2, yy1, yy2;

    if( m_usePin[moduleID1] )//false
    {
        xx1 = x[2 * moduleID1]     + m_pDB->m_pins[ pinID1 ].xOff;
        yy1 = x[2 * moduleID1 + 1] + m_pDB->m_pins[ pinID1 ].yOff;
    } else {
        xx1 = x[2 * moduleID1];
        yy1 = x[2 * moduleID1 + 1];
    }

    if ( m_usePin[moduleID2] )
    {
        xx2 = x[2 * moduleID2]     + m_pDB->m_pins[ pinID2 ].xOff;
        yy2 = x[2 * moduleID2 + 1] + m_pDB->m_pins[ pinID2 ].yOff;
    } else {
        xx2 = x[2 * moduleID2];
        yy2 = x[2 * moduleID2 + 1];
    }

    if ( m_pDB->m_modules[moduleID1].m_isFixed == false )
    {
        sg_getFirstModuleSubWireGrad(moduleID1, moduleID2, xx1, xx2, yy1, yy2,isPart);
    }

    if ( m_pDB->m_modules[moduleID2].m_isFixed == false )
    {
        sg_getFirstModuleSubWireGrad(moduleID2, moduleID1, xx2, xx1, yy2, yy1,isPart);
    }
}

void MyNLP::sg_updateSubGradWire_netBig_s(const int &netID, const double *x,bool isPart )
{
    double maxX = -10000;
    double maxY = -10000;
    double minX = DBL_MAX;
    double minY = DBL_MAX;

    int maxIDx, maxIDy, minIDx, minIDy;
    int pinID, moduleID;
    // To find the max and min area of the net
    for (unsigned int j = 0; j < m_pDB->m_nets[netID].size(); j++)
    {
        pinID = m_pDB->m_nets[netID][j];
        moduleID = m_pDB->m_pins[pinID].moduleId;

        double moduleX = x[2 * moduleID];
        double moduleY = x[2 * moduleID + 1];

        if (m_usePin[pinID] == true)
        {
            moduleX = x[2 * moduleID]     + m_pDB->m_pins[pinID].xOff;
            moduleY = x[2 * moduleID + 1] + m_pDB->m_pins[pinID].yOff;
        }

        if (moduleX >= maxX)
        {
            maxX = moduleX;
            maxIDx = moduleID;
        }
        if (moduleX <= minX)
        {
            minX = moduleX;
            minIDx = moduleID;
        }
        if (moduleY >= maxY)
        {
            maxY = moduleY;
            maxIDy = moduleID;
        }
        if (moduleY <= minY)
        {
            minY = moduleY;
            minIDy = moduleID;
        }
    }
    // To calulate the subgradient of the module
    //  int netSize = m_pDB->m_nets[netID].size();

    if (m_pDB->m_modules[maxIDx].m_isFixed == false)
    {
        sg_getXModuleSubWireGrad(maxIDx, minIDx, maxX, minX,isPart);
        //  sub_grad_wire[2 * maxIDx] += 1;
    }

    if (m_pDB->m_modules[minIDx].m_isFixed == false)
    {
        sg_getXModuleSubWireGrad(minIDx, maxIDx, minX, maxX,isPart);
        // sub_grad_wire[2 * minIDx] -= 1;
    }

    if (m_pDB->m_modules[maxIDy].m_isFixed == false)
    {
        sg_getYModuleSubWireGrad(maxIDy, minIDy, maxY, minY,isPart);
        // sub_grad_wire[2 * maxIDy + 1] += 1;
    }

    if (m_pDB->m_modules[minIDy].m_isFixed == false)
    {
        sg_getYModuleSubWireGrad(minIDy, maxIDy, minY, maxY,isPart);
        // sub_grad_wire[2 * minIDy + 1] -= 1;
    }



}

void MyNLP::sg_getFirstModuleSubWireGrad(const int& id1, const int& id2, const double& x1, const double& x2, const double& y1, const double& y2,bool isPart)
{

    sg_getXModuleSubWireGrad(id1, id2, x1, x2,isPart);
    sg_getYModuleSubWireGrad(id1, id2, y1, y2,isPart);
}

void MyNLP::sg_getXModuleSubWireGrad(const int& id1, const int& id2, const double& x1, const double& x2,bool isPart )
{
    if(isPart == false){
        if ( fabs( x1 - x2 ) < sgEps) {

            if (m_pDB->m_modules[id1].m_cx > m_pDB->m_modules[id2].m_cx )
            {
                sub_grad_wire[2 * id1] += getRandNum(0.5, 1.0);
                //sub_grad_wire[2 * id1] -= 0.0;
            } else if (m_pDB->m_modules[id1].m_cx < m_pDB->m_modules[id2].m_cx ) {
                sub_grad_wire[2 * id1] -= getRandNum(0.5, 1.0);
                //sub_grad_wire[2 * id1] += 0.0;
            } else {
                // TODO: compare the area to determine the direction ??
                //  if (m_pDB->m_modules[id1].m_area < m_pDB->m_modules[id2].m_area)
                {
                    sub_grad_wire[2 * id1] += 1;
                }
            }

        } else if (x1 > x2) {
            sub_grad_wire[2 * id1] += 1;
        } else {
            sub_grad_wire[2 * id1] -= 1;
        }
    }
    else{
        if ( fabs( x1 - x2 ) < sgEps) {

            if (m_pDB->m_modules[id1].m_cx > m_pDB->m_modules[id2].m_cx )
            {
                part_grad_wire[2 * id1] += getRandNum(0.5, 1.0);
                //sub_grad_wire[2 * id1] -= 0.0;
            } else if (m_pDB->m_modules[id1].m_cx < m_pDB->m_modules[id2].m_cx ) {
                part_grad_wire[2 * id1] -= getRandNum(0.5, 1.0);
                //sub_grad_wire[2 * id1] += 0.0;
            } else {
                // TODO: compare the area to determine the direction ??
                //  if (m_pDB->m_modules[id1].m_area < m_pDB->m_modules[id2].m_area)
                {
                    part_grad_wire[2 * id1] += 1;
                }
            }

        } else if (x1 > x2) {
            part_grad_wire[2 * id1] += 1;
        } else {
            part_grad_wire[2 * id1] -= 1;
        }

    }
}

void MyNLP::sg_getYModuleSubWireGrad(const int& id1, const int& id2, const double& y1, const double& y2,bool isPart)
{
    // for the precision
    if(isPart == false){
        if ( fabs( y1 - y2 ) < sgEps ) {
            if (m_pDB->m_modules[id1].m_cy > m_pDB->m_modules[id2].m_cy ) {
                sub_grad_wire[2 * id1 + 1] += getRandNum(0.5, 1.0);
                //sub_grad_wire[2 * id1 + 1] -= 0.0;
            } else if (m_pDB->m_modules[id1].m_cy < m_pDB->m_modules[id2].m_cy) {
                sub_grad_wire[2 * id1 + 1] -= getRandNum(0.5, 1.0);
                //sub_grad_wire[2 * id1 + 1] += 0.0;
            } else {
                //TODO: compare the area to determine the direction
                //  if (m_pDB->m_modules[id1].m_area < m_pDB->m_modules[id2].m_area)
                {
                    sub_grad_wire[2 * id1 + 1] += 1;
                }
            }
        } else if (y1 > y2) {
            sub_grad_wire[2 * id1 + 1] += 1;
        } else {
            sub_grad_wire[2 * id1 + 1] -= 1;
        }
    }
    else{
        if ( fabs( y1 - y2 ) < sgEps ) {
            if (m_pDB->m_modules[id1].m_cy > m_pDB->m_modules[id2].m_cy ) {
                part_grad_wire[2 * id1 + 1] += getRandNum(0.5, 1.0);
                //sub_grad_wire[2 * id1 + 1] -= 0.0;
            } else if (m_pDB->m_modules[id1].m_cy < m_pDB->m_modules[id2].m_cy) {
                part_grad_wire[2 * id1 + 1] -= getRandNum(0.5, 1.0);
                //sub_grad_wire[2 * id1 + 1] += 0.0;
            } else {
                //TODO: compare the area to determine the direction
                //  if (m_pDB->m_modules[id1].m_area < m_pDB->m_modules[id2].m_area)
                {
                    part_grad_wire[2 * id1 + 1] += 1;
                }
            }
        } else if (y1 > y2) {
            part_grad_wire[2 * id1 + 1] += 1;
        } else {
            part_grad_wire[2 * id1 + 1] -= 1;
        }
    }
}

/**********************************************flpeng added***********************************************************/

/********************************************hdgao added***************************************************************/
bool   MyNLP::SLSolve_nets_cluster(double wWire, double target_density, int currentLevel)
{

    //  test(); // added by hdgao

    double givenTargetUtil = m_targetUtil;    // m_targetUtil = -1
    m_currentStep = param.step;

    m_targetUtil += 0.05;
    if (m_targetUtil > 1.0)
        m_targetUtil = 1.0;

    double time_start = seconds();
    char filename[100];    // for gnuplot

    int n = 2 * m_pDB->m_modules.size();  //for the total num of modules

    // calculate the goal of the utilization
    double designUtil = m_pDB->m_totalModuleArea / m_pDB->m_totalFreeSpace;

    if (param.bShow)
        printf( "hdgao INFO: Design utilization: %f\n", designUtil);

    if (m_targetUtil > 0)
    {  // has give a utilization
        double lowest = designUtil + 0.05;
        if (m_targetUtil < lowest)
        {
            if (param.bShow)
            {
                printf("hdgao WARNING: Target utilization (%f) is too low\n", m_targetUtil);
                printf("hdgao          Set target utilization to %f \n", lowest);
            }
            m_targetUtil = lowest;
        }
    } else { // no given utilization
        printf("hdgao WARNING: No given target utilization. Distribute blocks evenly.\n");
        m_targetUtil = designUtil + 0.05;
        if (m_targetUtil > 1.0)
            m_targetUtil = 1.0;
    }

    if (param.bShow)
        printf("hdgao DBIN: Target utilization: %f\n", m_targetUtil);



    double* lastPosition = new double[n];
    memset(lastPosition, 0, sizeof(double) * n);


    //calculate the subgrad of f and the direction and allocate memeory
    double* sub_grad_f = new double [ n];
    double* last_sub_grad_f = new double [ n ]; // for computing CG-direction;
    double* real_last_grad_f = new double [ n];
    double* old_grad_f = new double[n];
    memset( sub_grad_f, 0, sizeof(double)*n );
    memset( last_sub_grad_f, 0, sizeof(double)*n );
    memset(real_last_grad_f, 0, sizeof(double) * n);
    memset(old_grad_f, 0, sizeof(double) *n);



    int densityGridSize = 10;
    double objValue;

    // init
    CreatePotentialGrid();
    CreateDensityGrid(densityGridSize);
    // fixed value
    UpdateDensityGridSpace(n, x);// ? @hdgao2504
    UpdatePotentialGridBase(x);
    UpdateExpBinPotential(m_targetUtil);
    // update when modules move
    UpdatePotentialGrid(x);
    UpdateDensityGrid(n, x);

    _weightWire = 1.0;

    density = GetDensityPanelty();// to evaluate overlap area  @2704hdgao
    sg_eval_f(x, objValue);
    sg_eval_grad_f(n, x, sub_grad_f);
    AdjustForce(n, x, sub_grad_wire, sub_grad_potential); // ? modify sub_grad_wire and sub_grad_potential

    double totalSubWireGradient = 0.0;
    double totalSubPotentialGradient = 0.0;
    for (int i = 0; i < n; i++) {
        totalSubWireGradient      += fabs(sub_grad_wire[i]);
        totalSubPotentialGradient += fabs(sub_grad_potential[i]);
    }

    _weightDensity = 1.0*1.02;// hdgao changed
    _weightWire = wWire * totalSubPotentialGradient / totalSubWireGradient; //* param.sg_namdaIncFactor;


    double   beta;   // for the CG method
    double nnbReal = GetNonZeroDensityGridPercent();
    double maxDen = GetMaxDensity();
    double totalOverDen = GetTotalOverDensity();
    double totalOverDenLB = GetTotalOverDensityLB();
    double totalOverPotential = GetTotalOverPotential();

    if (param.bShow)
    {
        printf(" %d-%2d HPWL= %.0f\tDen= %.2f %.2f %.2f %.2f NNB= %.2f Dcost= %4.1f%%  WireW= %.0f",
               currentLevel, m_ite, m_pDB->CalcHPWL(), maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
               nnbReal, density * _weightDensity / objValue * 100.0, _weightWire
               );
    } else {
        printf( " %d-%2d HPWL= %.0f \t", currentLevel, m_ite, m_pDB->CalcHPWL() );
    }
    fflush (stdout);

    if( param.bShow )
    {
        sprintf( filename, "SLfig%d-%d.plt", currentLevel, m_ite );
        m_pDB->OutputGnuplotFigure( filename, false, false );
    }


    double lastTotalOver = 0.0;
    double lastTotalOverPotential = DBL_MAX;
    double over = totalOverDen;
    //    int    totalIte = 0;
    bool   hasBestLegalSol = false;
    double bestLegalWL = DBL_MAX;
    int    lookAheadLegalCount = 0;
    double totalLegalTime = 0.0;
    bool   startDecreasing = false;
    int checkStep = 5;
    int LALnoGoodCount = 0;

    srand(time(NULL));
    double ssModify = 1.0;
    double stepSize = 1; // 1 - > 2
    double* x_bak = new double[n];


    //    nets_cluster_analysics();
    NCA_deeper();


    nets_cluster_modules_bounds.resize(nets_cluster.size());
    for(unsigned int i = 0;i < nets_cluster.size();i++){
        getBlackModule(i);
    }



    printf("\n======================================================================================\n");
    printf("cg begin\n");
    int      maxIte = int( log(_weightWire)/log(m_incFactor) );   // max iterator
    bool is_2D = false;
    for(int ite = 0;ite < maxIte;ite++){
        //        while(true){
        //            int update_counts = 0;
        //        double OneiteTime = seconds();
        //            m_currentStep = param.sg_stepsize;//param.sg_stepsize;  param.step
        m_currentStep = 0.3;
        vector<double> SL_grad_f;
        SL_grad_f.resize(n,0.0);


        for(int i = 0; i < n;i++)
            x_bak[i] = x[i];


        double last_objValue = DBL_MAX;
        double now_objValue = DBL_MAX;
        double old_obj,last_obj_value,obj_value;

        unsigned int nets_cluster_size = nets_cluster.size();
        double obj_before_OneIte,obj_after_OneIte;
        sg_eval_f(x,obj_before_OneIte);
        //            printf("obj_before_OneIte = %f\n",obj_after_OneIte);
        //*******************************************************CG begin***********************************************************************//
        for(unsigned int i_nets_cluster = 0; i_nets_cluster < nets_cluster_size;i_nets_cluster++){
            double OneClusterTime = seconds();
            vector<double> last_grad_f,grad_f;
            last_grad_f.resize(n,0.0);
            grad_f.resize(n,0.0);
            vector<int> RegionModule;

            updateBlackBounds(i_nets_cluster);
            part_eval_grad_f(x,i_nets_cluster,RegionModule,grad_f); // get grad

            if(is_2D)
                PartAdjustForce_2D(grad_f,RegionModule);
            else
                PartAdjustForce(grad_f,RegionModule);    // m_gridDensity is  irrelevant to sub_grad_f

            for(int i_cg = 0;true;i_cg++){// 5 -> 3

                last_grad_f.assign(grad_f.begin(),grad_f.end());
                updateBlackBounds(i_nets_cluster);
                part_eval_grad_f(x,i_nets_cluster,RegionModule,grad_f,false); // get grad
                //                   part_eval_grad_f(x,i_nets_cluster,RegionModule,grad_f);
                if(is_2D)
                    PartAdjustForce_2D(grad_f,RegionModule);
                else
                    PartAdjustForce(grad_f,RegionModule);    // m_gridDensity is  irrelevant to sub_grad_f

                if(i_cg % checkStep == 0){
                    old_obj = last_obj_value;
                    sg_eval_f(x,last_obj_value);
                }
                if(i_cg % checkStep == 0){
                    printf(".");
                    if(i_cg % (2 * checkStep) == 0 ){
                        UpdateBlockPosition(x);
                        if( sg_CalcHPWL(x) > bestLegalWL ){//??  outside
                            printf("get best legalWL\n");
                            break;
                        }
                    }

                    UpdateDensityGrid( n, x );  // find the exact bin density
                    totalOverDen = GetTotalOverDensity();
                    totalOverDenLB = GetTotalOverDensityLB();
                    totalOverPotential = GetTotalOverPotential();

                    double lastTotalOver = over;
                    double over = totalOverPotential;//min( totalOverPotential, totalOverDen ); // TEST  min

                    if( !startDecreasing    // outside
                            && over < lastTotalOver
                            && i_cg >= 5  // hyperparameter
                            && ite >= 1){  // outside
                        printf( ">>" );
                        fflush( stdout );
                        startDecreasing = true;
                    }
                    if( startDecreasing && over < target_density ) // ???
                        break;
                    if( obj_value >= param.precision * old_obj){
                        break;
                    }


                }
                //calculate d_k
                if (i_cg == 0){// first
                    for (int i = 0; i < RegionModule.size(); i++){
                        int moduleId = RegionModule[i];
                        grad_f[2*moduleId]   = -grad_f[2*moduleId];
                        grad_f[2*moduleId+1] = -grad_f[2*moduleId+1];
                    }
                }
                else {
                    // conjugate sub gradient direction
                    // FindBeta(n, sub_grad_f, last_sub_grad_f, beta);
                    double beta1,beta2,beta;
                    if(is_2D)
                        PartFindBeta_2D(grad_f, last_grad_f,RegionModule, beta1,beta2 );
                    else{
                        PartFindBeta(grad_f, last_grad_f,RegionModule, beta );
                        beta1 = beta2 = beta;
                    }

                    for (int i = 0; i < RegionModule.size(); i++){
                        int moduleId = RegionModule[i];
                        grad_f[2*moduleId]   = -grad_f[2*moduleId]   + beta1 * last_grad_f[2*moduleId];
                        grad_f[2*moduleId+1] = -grad_f[2*moduleId+1] + beta2 * last_grad_f[2*moduleId+1];

                    }
                }

                //update x , the step size
                double stepSize1,stepSize2,stepSize;
                if(is_2D)
                    PartLineSearch_2D(grad_f,RegionModule,stepSize1,stepSize2);
                else{
                    PartLineSearch(grad_f,RegionModule,stepSize);
                    stepSize1 = stepSize2 = stepSize;
                }

                for(int i = 0;i < RegionModule.size();i++){
                    int moduleId = RegionModule[i];
                    x[2*moduleId]   += grad_f[2*moduleId  ] * stepSize1;
                    x[2*moduleId+1] += grad_f[2*moduleId+1] * stepSize2;
                    SL_grad_f[2*moduleId]   += grad_f[2*moduleId];
                    SL_grad_f[2*moduleId+1] += grad_f[2*moduleId+1];
                }

                int n = 2*m_pDB->m_modules.size();
                BoundX(n, x, x_l, x_u);
                //       double now_obj;
                //       sg_eval_f(x,now_obj);
                //       if(now_obj < old_obj){
                //            printf("down\t");
                //            decrease_count++;
                //       }
                vector<int> right_mId;
                for(int i = 0; i < n/2;i++){
                    if( x_bak[2*i] != x[2*i] || x_bak[2*i + 1] != x[2*i + 1])
                        right_mId.push_back(i);
                }

                part_UpdatePotentialGrid(x,x_bak,right_mId,right_mId.size());

                //                 updateBlackBounds(i_nets_cluster);
                for(int i = 0; i < right_mId.size();i++){
                    x_bak[ 2*right_mId[i] ]   = x[ 2*right_mId[i] ];
                    x_bak[ 2*right_mId[i]+1 ] = x[ 2*right_mId[i]+1 ];
                }



            }// end optimize one group net

            UpdateBlockPosition(x);
            UpdateDensityGrid(n,x);

        }// end traversal net_cluster
        //*******************************************************CG   end***********************************************************************//

        //            printf("\ntraversal time = %f\n",s econds() - OneiteTime);
        sg_eval_f(x,obj_after_OneIte);
        if(obj_after_OneIte - obj_before_OneIte <= 0)
            printf(" down! down! down! down! down! down! down! down! down! down! down! down!\n");
        else
            printf("up! up! up! up! up! up! up! up! up! up! up! up! up! up! up! up! up! up! \n");
        //            getchar();
#if 0
        int merge_step = 1;
        if(ite % merge_step == 0){//
            double SL_step = 0;
            AdjustForce(n,x,SL_grad_f);
            LineSearch(n, x, SL_grad_f, SL_step);
            for(unsigned int i = 0; i < m_pDB->m_modules.size();i++){
                x[2*i]   += SL_step*SL_grad_f[2*i];
                x[2*i + 1] += SL_step*SL_grad_f[2*i + 1];
                lastPosition[2*i] = x[2*i];
                lastPosition[2*i+1] = x[2*i+1];
                assert(x[2*i] == x[2*i]);
                assert(x[2*i+1] == x[2*i+1]);
            }

            BoundX(n, x, x_l, x_u);
            UpdateBlockPosition(x);
            UpdateDensityGrid(n,x);
            UpdatePotentialGrid(x);
        }
#endif
        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();
        sg_eval_f(x, objValue);
        if(param.bShow){
            printf( "%d-%2d HPWL= %.0f\maxDen= %.2f OD= %.4f ODLB=%.4f OP=%.4f  LTime= %.1fm  objValue=%.0f, WireW= %.0f \n",
                    currentLevel, ite+1, m_pDB->CalcHPWL(),
                    maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                    double(seconds()-time_start)/60.0,
                    objValue,
                    _weightWire
                    );
        }

        if( param.bShow )
        {
            sprintf( filename, "SLfig%d-%d.plt", currentLevel, ite );
            m_pDB->OutputGnuplotFigure( filename, false, false );
        }

        if(ite % checkStep == 0 ){//record objValue per checkStep

            last_objValue = now_objValue;
            sg_eval_f(x, objValue);
            now_objValue = objValue;

            //                UpdateBlockPosition(x);// update modules positions per checkStep
            //                last_wirelength = now_wirelength;
            //                now_wirelength = m_pDB->CalcHPWL();

        }

        if( ite % checkStep == 0 )
        {
            //              printf( "." );
            fflush( stdout );

            if( ite % (2 * checkStep) == 0 )
            {
                UpdateBlockPosition( x );   // update to placeDB
                if( m_pDB->CalcHPWL() > bestLegalWL )   // gWL > LAL-WL
                {
                    printf( "BestLegalWL Break!\n" );
                    fflush( stdout );
                    break;
                }
            }

            UpdateDensityGrid( n, x );  // find the exact bin density
            totalOverDen = GetTotalOverDensity();
            totalOverDenLB = GetTotalOverDensityLB();
            totalOverPotential = GetTotalOverPotential();

            lastTotalOver = over;
            over = min( totalOverPotential, totalOverDen ); // TEST

            if( !startDecreasing
                    && over < lastTotalOver
                    && ite >= 6 )
            {
                printf( ">>" );
                fflush( stdout );
                startDecreasing = true;
            }

            // 2005-03-11: meet the constraint
            if( startDecreasing && over < target_density ){
                printf("meet target_density %f",target_density);
                break;
            }

            // Cannot further improve
            if( now_objValue >= param.precision * last_objValue )
            {
                printf("Cannot further improve = _ =\n");
                for (int i = 0; i < n; i++) {
                    x[i] = lastPosition[i];
                }
                break;
            }


        } // end of ite == checkStep

        if(ite%1 == 0) {

            for (int i = 0; i < n; i++) {
                lastPosition[i] = x[i];
            }
        }


        //        }// end while true

        UpdateDensityGrid(n,x);
        UpdateBlockPosition(x);

        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();

        if (ite >= 1 && m_lookAheadLegalization && over < target_density + 0.1 )
        {
            UpdateBlockPosition(x);
            double hpwl = m_pDB->CalcHPWL();
            if (hpwl > bestLegalWL)
            {
                printf("Stop. Good enough.\n");
                break;
            }

            lookAheadLegalCount++;
            double oldWL = hpwl;
            CTetrisLegal legal(*m_pDB);

            double scale = 0.85;
            if (givenTargetUtil < 1.0 && givenTargetUtil > 0)
            {
                scale = 0.9;
            }

            double legalStart = seconds();
            bool   bLegal = legal.Solve( givenTargetUtil, false, false, scale);
            double legalTime = seconds() - legalStart;

            totalLegalTime += legalTime;

            if (param.bShow)
            {
                printf("LAL Time: %.2f\n", legalTime);
            }

            if (bLegal)
            {
                double WL = m_pDB->GetHPWLdensity(givenTargetUtil);
                if (param.bShow)
                {
                    m_pDB->ShowDensityInfo();
                }

                if (WL < bestLegalWL)
                {
                    LALnoGoodCount = 0;
                    if (param.bShow)
                    {
                        printf("SAVE BEST! (HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n",
                               m_pDB->GetHPWLp2p(), WL, (WL - oldWL) / oldWL * 100);
                    }
                    bestLegalWL = WL;
                    hasBestLegalSol = true;

                    for (int i = 0; i < (int) m_pDB->m_modules.size(); i++)
                    {
                        xBest[2 * i]     = m_pDB->m_modules[i].m_cx;
                        xBest[2 * i + 1] = m_pDB->m_modules[i].m_cy;
                        //                        x[2 * i]     = m_pDB->m_modules[i].m_cx;
                        //                        x[2 * i + 1] = m_pDB->m_modules[i].m_cy;
                    }
                } else {
                    if( param.bShow )
                        printf( "(HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n", m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    if( (WL-oldWL)/oldWL < 0.01 )
                    {
                        if( param.bShow )
                            printf( "Stop. Good enough\n" );
                        break;
                    }
                    LALnoGoodCount++;
                    if( LALnoGoodCount >= 2 )
                        break;
                }
            }

        }
        //end of m_lookAheadLegalization
        if (ite >= 1)
        {
            if ( m_lookAheadLegalization && startDecreasing && over < target_density + 0.01)
            {
                printf("Meet constraint!\n");
                break;
            }
            if (!m_lookAheadLegalization && startDecreasing && over < target_density + 0.03)
            {
                printf("Meet constraint!\n");
                break;
            }

            if (ite > 3 && totalOverPotential > lastTotalOverPotential && totalOverPotential < 1.4)
            {
                printf("Cannot further reduce!\n");
                break;
            }
        }
        // reduce _weightwire instead of _weightDensity  @hdgao
        /*
        if (over > 1.5)
        {
            _weightWire /= 1.5;
        } else if (over > 1.1) {
            _weightWire /= 1.5;
        } else {
            _weightWire /= 1.3;
        }
*/
        //       _weightDensity *=2;  //  _weightWire /=2;
        _weightWire /= m_incFactor;
        lastTotalOverPotential = totalOverPotential;
        //        printf("TOTAL one ite %d time = %f\n",ite,seconds() - OneiteTime);
    }// end of "for (ite = 0; ite < maxIte; ite++)"

    if (hasBestLegalSol)
    {
        memcpy(x, xBest, sizeof(double) * n);
    }

    UpdateBlockPosition(x);

    if( lookAheadLegalCount > 0 && param.bShow )
    {
        printf( "LAL: Total Count: %d\n", lookAheadLegalCount );
        printf( "LAL: Total CPU: %.2f\n", totalLegalTime );
        sprintf( filename, "util-global.dat" );
        CPlaceBin placeBin( *m_pDB );
        placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
        placeBin.OutputBinUtil( filename );
    }

    delete [] lastPosition;

    delete [] sub_grad_f;
    delete [] last_sub_grad_f;
    delete [] real_last_grad_f;
    delete [] old_grad_f;
    return hasBestLegalSol;
}

bool   MyNLP::SLSolve_bak(double wWire, double target_density, int currentLevel)
{



    double givenTargetUtil = m_targetUtil;    // m_targetUtil = -1
    m_currentStep = param.step;

    m_targetUtil += 0.05;
    if (m_targetUtil > 1.0)
        m_targetUtil = 1.0;

    double time_start = seconds();
    char filename[100];    // for gnuplot

    int n = 2 * m_pDB->m_modules.size();  //for the total num of modules

    // calculate the goal of the utilization
    double designUtil = m_pDB->m_totalModuleArea / m_pDB->m_totalFreeSpace;

    if (param.bShow)
        printf( "hdgao INFO: Design utilization: %f\n", designUtil);

    if (m_targetUtil > 0)
    {  // has give a utilization
        double lowest = designUtil + 0.05;
        if (m_targetUtil < lowest)
        {
            if (param.bShow)
            {
                printf("hdgao WARNING: Target utilization (%f) is too low\n", m_targetUtil);
                printf("hdgao          Set target utilization to %f \n", lowest);
            }
            m_targetUtil = lowest;
        }
    } else { // no given utilization
        printf("flpeng WARNING: No given target utilization. Distribute blocks evenly.\n");
        m_targetUtil = designUtil + 0.05;
        if (m_targetUtil > 1.0)
            m_targetUtil = 1.0;
    }

    if (param.bShow)
        printf("flpeng DBIN: Target utilization: %f\n", m_targetUtil);



    double* lastPosition = new double[n];
    memset(lastPosition, 0, sizeof(double) * n);


    //calculate the subgrad of f and the direction and allocate memeory
    double* sub_grad_f = new double [ n];
    double* last_sub_grad_f = new double [ n ]; // for computing CG-direction;
    double* real_last_grad_f = new double [ n];
    double* old_grad_f = new double[n];
    memset( sub_grad_f, 0, sizeof(double)*n );
    memset( last_sub_grad_f, 0, sizeof(double)*n );
    memset(real_last_grad_f, 0, sizeof(double) * n);
    memset(old_grad_f, 0, sizeof(double) *n);

    double* SL_grad_f = new double [n];
    memset( SL_grad_f, 0, sizeof(double)*n );


    int densityGridSize = 10;
    double objValue;

    // init
    CreatePotentialGrid();
    CreateDensityGrid(densityGridSize);
    // fixed value
    UpdateDensityGridSpace(n, x);// ? @hdgao2504
    UpdatePotentialGridBase(x);
    UpdateExpBinPotential(m_targetUtil);
    // update when modules move
    UpdatePotentialGrid(x);
    UpdateDensityGrid(n, x);

    _weightWire = 1.0;

    density = GetDensityPanelty();// to evaluate overlap area  @2704hdgao
    sg_eval_f(x, objValue);
    sg_eval_grad_f(n, x, sub_grad_f);
    AdjustForce(n, x, sub_grad_wire, sub_grad_potential); // ? modify sub_grad_wire and sub_grad_potential

    double totalSubWireGradient = 0.0;
    double totalSubPotentialGradient = 0.0;
    for (int i = 0; i < n; i++) {
        totalSubWireGradient      += fabs(sub_grad_wire[i]);
        totalSubPotentialGradient += fabs(sub_grad_potential[i]);
    }

    _weightDensity = 1.0;
    _weightWire = wWire * totalSubPotentialGradient / totalSubWireGradient * param.sg_namdaIncFactor;// * param.sg_namdaIncFactor

    int      maxIte = 50;   // max iterator
    double   beta;   // for the CG method
    double nnbReal = GetNonZeroDensityGridPercent();
    double maxDen = GetMaxDensity();
    double totalOverDen = GetTotalOverDensity();
    double totalOverDenLB = GetTotalOverDensityLB();
    double totalOverPotential = GetTotalOverPotential();

    if (param.bShow)
    {
        printf(" %d-%2d HPWL= %.0f\tDen= %.2f %.2f %.2f %.2f NNB= %.2f Dcost= %4.1f%%  WireW= %.0f",
               currentLevel, m_ite, m_pDB->CalcHPWL(), maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
               nnbReal, density * _weightDensity / objValue * 100.0, _weightWire
               );
    } else {
        printf( " %d-%2d HPWL= %.0f \t", currentLevel, m_ite, m_pDB->CalcHPWL() );
    }
    fflush (stdout);

    if( param.bShow )
    {
        sprintf( filename, "SLfig%d-%d.plt", currentLevel, m_ite );
        m_pDB->OutputGnuplotFigure( filename, false, false );
    }


    double lastTotalOver = 0.0;
    double lastTotalOverPotential = DBL_MAX;
    double over = totalOverDen;
    //    int    totalIte = 0;
    bool   hasBestLegalSol = false;
    double bestLegalWL = DBL_MAX;
    int    lookAheadLegalCount = 0;
    double totalLegalTime = 0.0;
    bool   startDecreasing = false;
    int checkStep = 1;
    int LALnoGoodCount = 0;

    srand(time(NULL));
    double ssModify = 1.0;
    double stepSize = 1; // 1 - > 2
    double* x_bak = new double[n];






    printf("\n======================================================================================\n");
    printf("cg begin\n");
    for(int ite = 0;ite < maxIte;ite++){
        //        while(true){
        //            int update_counts = 0;
        m_currentStep = param.sg_stepsize;
        memset( SL_grad_f, 0, sizeof(double)*n );

        for(int i = 0; i < n;i++)
            x_bak[i] = x[i];


        double last_objValue = DBL_MAX;
        double now_objValue = DBL_MAX;

        //            double last_wirelength = DBL_MAX;
        //           double now_wirelength  = DBL_MAX;
        //            double best_wirelenth = DBL_MAX;

        //  update modules positions
        int num_nets = 60;
        bool nets_flag = true;
        unsigned int total_nets = m_pDB->m_nets.size();
        for(unsigned int i_nets = 0; nets_flag ;i_nets += num_nets){
            if(i_nets + num_nets -1 >= total_nets){
                num_nets = total_nets - i_nets;
                nets_flag = false;
            }
            int current_module_size = 0;
            for(int i = 0;i<num_nets;i++){
                current_module_size += m_pDB->m_nets[ i+i_nets ].size();
            }

            memset( sub_grad_f, 0, sizeof(double)*2*current_module_size);
            memset( last_sub_grad_f, 0, sizeof(double)*2*current_module_size);
            memset(real_last_grad_f, 0, sizeof(double) * 2 * current_module_size);
            memset(old_grad_f, 0, sizeof(double) * 2 * current_module_size);

            //                part_eval_f(x, objValue,i_nets);// get objValue
            part_eval_grad_f(x, sub_grad_f,i_nets,num_nets); // get grad
            AdjustForce(2*current_module_size,x,sub_grad_f);    // m_gridDensity is  irrelevant to sub_grad_f

            for(int i_cg = 0;i_cg < 5;i_cg++){// 5 -> 3
                swap(real_last_grad_f, old_grad_f);
                swap(last_sub_grad_f, sub_grad_f);
                part_eval_grad_f(x, sub_grad_f,i_nets,num_nets);
                AdjustForce(2*current_module_size,x,sub_grad_f);
                //calculate d_k
                if (i_cg == 0){// first
                    for (int i = 0; i < 2*current_module_size; i++){
                        old_grad_f[i] = sub_grad_f[i];
                        sub_grad_f[i] = -sub_grad_f[i];
                    }
                }
                else {
                    // conjugate sub gradient direction
                    sg_findBeta(2*current_module_size, sub_grad_f, last_sub_grad_f, real_last_grad_f, beta);
                    // FindBeta(n, sub_grad_f, last_sub_grad_f, beta);
                    for (int i = 0; i < 2*current_module_size; i++)
                    {
                        old_grad_f[i] = sub_grad_f[i];
                        sub_grad_f[i] = -sub_grad_f[i] + beta * last_sub_grad_f[i];
                    }
                }

                //calculate a_k , the step size

                LineSearch(2*current_module_size, x, sub_grad_f, stepSize);

                for(int j = 0,k = 0;j < num_nets;j++){
                    for(unsigned int i = 0;i < m_pDB->m_nets[i_nets + j].size();i++){
                        int pinId =  m_pDB->m_nets[i_nets + j][i];
                        int moduleId =  m_pDB->m_pins[pinId].moduleId;
                        x[2*moduleId]   += sub_grad_f[2*k]   * stepSize;
                        x[2*moduleId+1] += sub_grad_f[2*k+1] * stepSize;
                        // record the final part_grad
                        //  if(i_cg == 4){// record the last grad
                        SL_grad_f[2*moduleId]     += sub_grad_f[2*k];
                        SL_grad_f[2*moduleId+1]   += sub_grad_f[2*k+1];
                        k++;
                        //}
                    }
                }
                BoundX(n, x, x_l, x_u);

                vector<int> right_mId;
                for(int i = 0; i < n/2;i++){
                    if( x_bak[2*i] != x[2*i] || x_bak[2*i + 1] != x[2*i + 1])
                        right_mId.push_back(i);
                }
                if(0){
                    //                       part_UpdatePotentialGrid(x,x_bak,adjust_moduleId,current_module_size);
                    part_UpdatePotentialGrid(x,x_bak,right_mId,right_mId.size());
                    //                       part_UpdatePotentialGrid(x);
                    vector< vector<double> > m_gridPotential_bak(m_gridPotential);
                    for(int i = 0;i < m_potentialGridSize;i++){
                        assert(m_gridPotential_bak == m_gridPotential);
                    }
                    //                     m_gridPotential_bak.assign(m_gridPotential.begin(),m_gridPotential.end());
                    UpdatePotentialGrid(x);
                    int count = 0;
                    for(int i = 0;i < m_potentialGridSize;i++){
                        for(int j = 0;j < m_potentialGridSize;j++){
                            if( fabs( m_gridPotential_bak[i][j] - m_gridPotential[i][j] ) > 1e-5 ){
                                cout << i << "--" << j << "  " <<fabs( m_gridPotential_bak[i][j] - m_gridPotential[i][j] )<< endl;
                                count++;
                            }
                        }
                    }
                    if(count) cout <<"############" << count << " of " << current_module_size << endl;
                    //                       for(int i = 0; i < current_module_size;i++){
                    //                           int temp_mId = adjust_moduleId[i];
                    //                           x_bak[2*temp_mId] = x[2*temp_mId];
                    //                           x_bak[2*temp_mId + 1] = x[2*temp_mId + 1];
                    //                       }
                    for(int i = 0; i < right_mId.size();i++){
                        x_bak[ 2*right_mId[i] ]   = x[ 2*right_mId[i] ];
                        x_bak[ 2*right_mId[i]+1 ] = x[ 2*right_mId[i]+1 ];
                    }
                    if(count) getchar();
                    continue;
                }

                if(0){
                    printf("netID = %d;,net_size = %d:\n",i_nets,current_module_size);
                    for(int i = 0; i < current_module_size;i++){
                        printf("grad[%d] = %f,grad[%d] = %f;\t",2*i,sub_grad_f[2*i],
                                2*i+1,sub_grad_f[2*i+1]);
                    }
                    printf("\n");
                }

                part_UpdatePotentialGrid(x,x_bak,right_mId,right_mId.size());
                //                   UpdatePotentialGrid(x);    //time-consuming  0507

                for(int i = 0; i < right_mId.size();i++){
                    x_bak[ 2*right_mId[i] ]   = x[ 2*right_mId[i] ];
                    x_bak[ 2*right_mId[i]+1 ] = x[ 2*right_mId[i]+1 ];
                }
            }// end optimize one group net

            UpdateBlockPosition(x);
            UpdateDensityGrid(n,x);
            if(!param.bShow){
                maxDen = GetMaxDensity();
                totalOverDen = GetTotalOverDensity();
                totalOverDenLB = GetTotalOverDensityLB();
                totalOverPotential = GetTotalOverPotential();

                printf( " %d-%6d HPWL= %.0f\tDen= %.2f %.4f %.4f %.4f  LTime= %.1fm Dcost= %4.1f%% WireW= %.0f \n",
                        currentLevel, i_nets, m_pDB->CalcHPWL(),
                        maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                        double(seconds()-time_start)/60.0,
                        density*_weightDensity /objValue * 100.0,
                        _weightWire
                        );
            }

        }// end traversal nets

        int merge_step = 1;
        if(ite % merge_step == 0){//
            double SL_step = 0;
            AdjustForce(n,x,SL_grad_f);
            LineSearch(n, x, SL_grad_f, SL_step);
            for(unsigned int i = 0; i < m_pDB->m_modules.size();i++){
                x[2*i]   += SL_step*SL_grad_f[2*i];
                x[2*i + 1] += SL_step*SL_grad_f[2*i + 1];
                lastPosition[i] = x[i];
            }

            BoundX(n, x, x_l, x_u);
            UpdateBlockPosition(x);
            UpdateDensityGrid(n,x);
            UpdatePotentialGrid(x);
        }




        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();
        if(param.bShow){
            printf( " %d-%2d HPWL= %.0f\tDen= %.2f %.4f %.4f %.4f  LTime= %.1fm Dcost= %4.1f%% WireW= %.0f \n",
                    currentLevel, ite+1, m_pDB->CalcHPWL(),
                    maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                    double(seconds()-time_start)/60.0,
                    density*_weightDensity /objValue * 100.0,
                    _weightWire
                    );
        }
#if 1
        if( param.bShow )
        {
            sprintf( filename, "SLfig%d-%d.plt", currentLevel, ite );
            m_pDB->OutputGnuplotFigure( filename, false, false );
        }

#endif
        if(ite % checkStep == 0 ){//record objValue per checkStep
            last_objValue = now_objValue;
            sg_eval_f(x, objValue);
            now_objValue = objValue;

            //                UpdateBlockPosition(x);// update modules positions per checkStep
            //                last_wirelength = now_wirelength;
            //                now_wirelength = m_pDB->CalcHPWL();

        }
        if( ite % checkStep == 0 )
        {
            printf( "." );
            fflush( stdout );

            if( ite % (2 * checkStep) == 0 )
            {
                UpdateBlockPosition( x );   // update to placeDB
                if( m_pDB->CalcHPWL() > bestLegalWL )   // gWL > LAL-WL
                {
                    printf( "BestLegalWL Break!\n" );
                    fflush( stdout );
                    break;
                }
            }

            UpdateDensityGrid( n, x );  // find the exact bin density
            totalOverDen = GetTotalOverDensity();
            totalOverDenLB = GetTotalOverDensityLB();
            totalOverPotential = GetTotalOverPotential();

            lastTotalOver = over;
            over = min( totalOverPotential, totalOverDen ); // TEST

            if( !startDecreasing
                    && over < lastTotalOver
                    && ite >= 6 )
            {
                printf( ">>" );
                fflush( stdout );
                startDecreasing = true;
            }

            // 2005-03-11: meet the constraint
            if( startDecreasing && over < target_density )
                break;

            // Cannot further improve
            if( now_objValue >= param.precision * last_objValue )
            {
                printf("Cannot further improve = _ =\n");
                for (int i = 0; i < n; i++) {
                    x[i] = lastPosition[i];
                }
                break;
            }

        } // end of ite == checkStep
        if(ite%1 == 0) {
            for (int i = 0; i < n; i++) {
                lastPosition[i] = x[i];
            }
        }

        //        }// end while true

        UpdateDensityGrid(n,x);
        UpdateBlockPosition(x);

        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();

        if (ite >= 1 && m_lookAheadLegalization && over < target_density + 0.1 )
        {
            UpdateBlockPosition(x);
            double hpwl = m_pDB->CalcHPWL();
            if (hpwl > bestLegalWL)
            {
                printf("Stop. Good enough.\n");
                break;
            }

            lookAheadLegalCount++;
            double oldWL = hpwl;
            CTetrisLegal legal(*m_pDB);

            double scale = 0.85;
            if (givenTargetUtil < 1.0 && givenTargetUtil > 0)
            {
                scale = 0.9;
            }

            double legalStart = seconds();
            bool   bLegal = legal.Solve( givenTargetUtil, false, false, scale);
            double legalTime = seconds() - legalStart;

            totalLegalTime += legalTime;

            if (param.bShow)
            {
                printf("LAL Time: %.2f\n", legalTime);
            }

            if (bLegal)
            {
                double WL = m_pDB->GetHPWLdensity(givenTargetUtil);
                if (param.bShow)
                {
                    m_pDB->ShowDensityInfo();
                }

                if (WL < bestLegalWL)
                {
                    LALnoGoodCount = 0;
                    if (param.bShow)
                    {
                        printf("SAVE BEST! (HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n",
                               m_pDB->GetHPWLp2p(), WL, (WL - oldWL) / oldWL * 100);
                    }
                    bestLegalWL = WL;
                    hasBestLegalSol = true;

                    for (int i = 0; i < (int) m_pDB->m_modules.size(); i++)
                    {
                        xBest[2 * i]     = m_pDB->m_modules[i].m_cx;
                        xBest[2 * i + 1] = m_pDB->m_modules[i].m_cy;
                        //                        x[2 * i]     = m_pDB->m_modules[i].m_cx;
                        //                        x[2 * i + 1] = m_pDB->m_modules[i].m_cy;

                    }
                } else {
                    if( param.bShow )
                        printf( "(HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n", m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    if( (WL-oldWL)/oldWL < 0.01 )
                    {
                        if( param.bShow )
                            printf( "Stop. Good enough\n" );
                        break;
                    }
                    LALnoGoodCount++;
                    if( LALnoGoodCount >= 2 )
                        break;
                }
            }
        } //end of m_lookAheadLegalization
        if (ite >= 1)
        {
            if ( m_lookAheadLegalization && startDecreasing && over < target_density + 0.01)
            {
                printf("Meet constraint!\n");
                break;
            }
            if (!m_lookAheadLegalization && startDecreasing && over < target_density + 0.03)
            {
                printf("Meet constraint!\n");
                break;
            }

            if (ite > 3 && totalOverPotential > lastTotalOverPotential && totalOverPotential < 1.4)
            {
                printf("Cannot further reduce!\n");
                break;
            }
        }
        // reduce _weightwire instead of _weightDensity  @hdgao
        /*
        if (over > 1.5)
        {
            _weightWire /= 1.5;
        } else if (over > 1.1) {
            _weightWire /= 1.5;
        } else {
            _weightWire /= 1.3;
        }
*/
        _weightDensity *=2;  //  _weightWire /=2;

        lastTotalOverPotential = totalOverPotential;

    }// end of "for (ite = 0; ite < maxIte; ite++)"

    if (hasBestLegalSol)
    {
        memcpy(x, xBest, sizeof(double) * n);
    }

    UpdateBlockPosition(x);

    if( lookAheadLegalCount > 0 && param.bShow )
    {
        printf( "LAL: Total Count: %d\n", lookAheadLegalCount );
        printf( "LAL: Total CPU: %.2f\n", totalLegalTime );
        sprintf( filename, "util-global.dat" );
        CPlaceBin placeBin( *m_pDB );
        placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
        placeBin.OutputBinUtil( filename );
    }

    delete [] lastPosition;

    delete [] sub_grad_f;
    delete [] last_sub_grad_f;
    delete [] real_last_grad_f;
    delete [] old_grad_f;
    delete [] SL_grad_f;
    return hasBestLegalSol;
}

bool   MyNLP::SLSolve(double wWire, double target_density, int currentLevel)
{



    double givenTargetUtil = m_targetUtil;    // m_targetUtil = -1
    m_currentStep = param.step;

    m_targetUtil += 0.05;
    if (m_targetUtil > 1.0)
        m_targetUtil = 1.0;

    double time_start = seconds();
    char filename[100];    // for gnuplot

    int n = 2 * m_pDB->m_modules.size();  //for the total num of modules

    // calculate the goal of the utilization
    double designUtil = m_pDB->m_totalModuleArea / m_pDB->m_totalFreeSpace;

    if (param.bShow)
        printf( "hdgao INFO: Design utilization: %f\n", designUtil);

    if (m_targetUtil > 0)
    {  // has give a utilization
        double lowest = designUtil + 0.05;
        if (m_targetUtil < lowest)
        {
            if (param.bShow)
            {
                printf("hdgao WARNING: Target utilization (%f) is too low\n", m_targetUtil);
                printf("hdgao          Set target utilization to %f \n", lowest);
            }
            m_targetUtil = lowest;
        }
    } else { // no given utilization
        printf("flpeng WARNING: No given target utilization. Distribute blocks evenly.\n");
        m_targetUtil = designUtil + 0.05;
        if (m_targetUtil > 1.0)
            m_targetUtil = 1.0;
    }

    if (param.bShow)
        printf("flpeng DBIN: Target utilization: %f\n", m_targetUtil);



    double* lastPosition = new double[n];
    memset(lastPosition, 0, sizeof(double) * n);


    //calculate the subgrad of f and the direction and allocate memeory
    double* sub_grad_f = new double [ n];
    double* last_sub_grad_f = new double [ n ]; // for computing CG-direction;
    double* real_last_grad_f = new double [ n];
    double* old_grad_f = new double[n];
    memset( sub_grad_f, 0, sizeof(double)*n );
    memset( last_sub_grad_f, 0, sizeof(double)*n );
    memset(real_last_grad_f, 0, sizeof(double) * n);
    memset(old_grad_f, 0, sizeof(double) *n);

    double* SL_grad_f = new double [n];
    memset( SL_grad_f, 0, sizeof(double)*n );


    int densityGridSize = 10;
    double objValue;

    // init
    CreatePotentialGrid();
    CreateDensityGrid(densityGridSize);
    // fixed value
    UpdateDensityGridSpace(n, x);// ? @hdgao2504
    UpdatePotentialGridBase(x);
    UpdateExpBinPotential(m_targetUtil);
    // update when modules move
    UpdatePotentialGrid(x);
    UpdateDensityGrid(n, x);

    _weightWire = 1.0;

    density = GetDensityPanelty();// to evaluate overlap area  @2704hdgao
    sg_eval_f(x, objValue);
    sg_eval_grad_f(n, x, sub_grad_f);
    AdjustForce(n, x, sub_grad_wire, sub_grad_potential); // ? modify sub_grad_wire and sub_grad_potential

    double totalSubWireGradient = 0.0;
    double totalSubPotentialGradient = 0.0;
    for (int i = 0; i < n; i++) {
        totalSubWireGradient      += fabs(sub_grad_wire[i]);
        totalSubPotentialGradient += fabs(sub_grad_potential[i]);
    }

    _weightDensity = 1.0;
    _weightWire = wWire * totalSubPotentialGradient / totalSubWireGradient * param.sg_namdaIncFactor;// * param.sg_namdaIncFactor

    int      maxIte = 50;   // max iterator
    double   beta;   // for the CG method
    double nnbReal = GetNonZeroDensityGridPercent();
    double maxDen = GetMaxDensity();
    double totalOverDen = GetTotalOverDensity();
    double totalOverDenLB = GetTotalOverDensityLB();
    double totalOverPotential = GetTotalOverPotential();

    if (param.bShow)
    {
        printf(" %d-%2d HPWL= %.0f\tDen= %.2f %.2f %.2f %.2f NNB= %.2f Dcost= %4.1f%%  WireW= %.0f",
               currentLevel, m_ite, m_pDB->CalcHPWL(), maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
               nnbReal, density * _weightDensity / objValue * 100.0, _weightWire
               );
    } else {
        printf( " %d-%2d HPWL= %.0f \t", currentLevel, m_ite, m_pDB->CalcHPWL() );
    }
    fflush (stdout);

    if( param.bShow )
    {
        sprintf( filename, "SLfig%d-%d.plt", currentLevel, m_ite );
        m_pDB->OutputGnuplotFigure( filename, false, false );
    }


    double lastTotalOver = 0.0;
    double lastTotalOverPotential = DBL_MAX;
    double over = totalOverDen;
    //    int    totalIte = 0;
    bool   hasBestLegalSol = false;
    double bestLegalWL = DBL_MAX;
    int    lookAheadLegalCount = 0;
    double totalLegalTime = 0.0;
    bool   startDecreasing = false;
    int checkStep = 1;
    int LALnoGoodCount = 0;

    srand(time(NULL));
    double ssModify = 1.0;
    double stepSize = 1; // 1 - > 2
    double* x_bak = new double[n];






    printf("\n======================================================================================\n");
    printf("cg begin\n");
    for(int ite = 0;ite < maxIte;ite++){
        //        while(true){
        //            int update_counts = 0;
        m_currentStep = param.sg_stepsize;
        memset( SL_grad_f, 0, sizeof(double)*n );

        for(int i = 0; i < n;i++)
            x_bak[i] = x[i];


        double last_objValue = DBL_MAX;
        double now_objValue = DBL_MAX;

        //            double last_wirelength = DBL_MAX;
        //           double now_wirelength  = DBL_MAX;
        //            double best_wirelenth = DBL_MAX;

        //  update modules positions
        int num_nets = 60;
        bool nets_flag = true;
        unsigned int total_nets = m_pDB->m_nets.size();
        double obj1,obj2;
        sg_eval_f(x,obj1);
        for(unsigned int i_nets = 0; nets_flag ;i_nets += num_nets){
            if(i_nets + num_nets -1 >= total_nets){
                num_nets = total_nets - i_nets;
                nets_flag = false;
            }
            int current_module_size = 0;
            for(int i = 0;i<num_nets;i++){
                current_module_size += m_pDB->m_nets[ i+i_nets ].size();
            }

            memset( sub_grad_f, 0, sizeof(double)*2*current_module_size);
            memset( last_sub_grad_f, 0, sizeof(double)*2*current_module_size);
            memset(real_last_grad_f, 0, sizeof(double) * 2 * current_module_size);
            memset(old_grad_f, 0, sizeof(double) * 2 * current_module_size);

            //                part_eval_f(x, objValue,i_nets);// get objValue
            part_eval_grad_f(x, sub_grad_f,i_nets,num_nets); // get grad
            AdjustForce(2*current_module_size,x,sub_grad_f);    // m_gridDensity is  irrelevant to sub_grad_f

            for(int i_cg = 0;i_cg < 5;i_cg++){// 5 -> 3
                swap(real_last_grad_f, old_grad_f);
                swap(last_sub_grad_f, sub_grad_f);
                part_eval_grad_f(x, sub_grad_f,i_nets,num_nets);
                AdjustForce(2*current_module_size,x,sub_grad_f);
                //calculate d_k
                if (i_cg == 0){// first
                    for (int i = 0; i < 2*current_module_size; i++){
                        old_grad_f[i] = sub_grad_f[i];
                        sub_grad_f[i] = -sub_grad_f[i];
                    }
                }
                else {
                    // conjugate sub gradient direction
                    sg_findBeta(2*current_module_size, sub_grad_f, last_sub_grad_f, real_last_grad_f, beta);
                    // FindBeta(n, sub_grad_f, last_sub_grad_f, beta);
                    for (int i = 0; i < 2*current_module_size; i++)
                    {
                        old_grad_f[i] = sub_grad_f[i];
                        sub_grad_f[i] = -sub_grad_f[i] + beta * last_sub_grad_f[i];
                    }
                }

                //calculate a_k , the step size

                LineSearch(2*current_module_size, x, sub_grad_f, stepSize);

                for(int j = 0,k = 0;j < num_nets;j++){
                    for(unsigned int i = 0;i < m_pDB->m_nets[i_nets + j].size();i++){
                        int pinId =  m_pDB->m_nets[i_nets + j][i];
                        int moduleId =  m_pDB->m_pins[pinId].moduleId;
                        x[2*moduleId]   += sub_grad_f[2*k]   * stepSize;
                        x[2*moduleId+1] += sub_grad_f[2*k+1] * stepSize;
                        // record the final part_grad
                        //  if(i_cg == 4){// record the last grad
                        SL_grad_f[2*moduleId]     += sub_grad_f[2*k];
                        SL_grad_f[2*moduleId+1]   += sub_grad_f[2*k+1];
                        k++;
                        //}
                    }
                }
                BoundX(n, x, x_l, x_u);

                vector<int> right_mId;
                for(int i = 0; i < n/2;i++){
                    if( x_bak[2*i] != x[2*i] || x_bak[2*i + 1] != x[2*i + 1])
                        right_mId.push_back(i);
                }
                if(0){
                    //                       part_UpdatePotentialGrid(x,x_bak,adjust_moduleId,current_module_size);
                    part_UpdatePotentialGrid(x,x_bak,right_mId,right_mId.size());
                    //                       part_UpdatePotentialGrid(x);
                    vector< vector<double> > m_gridPotential_bak(m_gridPotential);
                    for(int i = 0;i < m_potentialGridSize;i++){
                        assert(m_gridPotential_bak == m_gridPotential);
                    }
                    //                     m_gridPotential_bak.assign(m_gridPotential.begin(),m_gridPotential.end());
                    UpdatePotentialGrid(x);
                    int count = 0;
                    for(int i = 0;i < m_potentialGridSize;i++){
                        for(int j = 0;j < m_potentialGridSize;j++){
                            if( fabs( m_gridPotential_bak[i][j] - m_gridPotential[i][j] ) > 1e-5 ){
                                cout << i << "--" << j << "  " <<fabs( m_gridPotential_bak[i][j] - m_gridPotential[i][j] )<< endl;
                                count++;
                            }
                        }
                    }
                    if(count) cout <<"############" << count << " of " << current_module_size << endl;
                    //                       for(int i = 0; i < current_module_size;i++){
                    //                           int temp_mId = adjust_moduleId[i];
                    //                           x_bak[2*temp_mId] = x[2*temp_mId];
                    //                           x_bak[2*temp_mId + 1] = x[2*temp_mId + 1];
                    //                       }
                    for(int i = 0; i < right_mId.size();i++){
                        x_bak[ 2*right_mId[i] ]   = x[ 2*right_mId[i] ];
                        x_bak[ 2*right_mId[i]+1 ] = x[ 2*right_mId[i]+1 ];
                    }
                    if(count) getchar();
                    continue;
                }

                if(0){
                    printf("netID = %d;,net_size = %d:\n",i_nets,current_module_size);
                    for(int i = 0; i < current_module_size;i++){
                        printf("grad[%d] = %f,grad[%d] = %f;\n",2*i,sub_grad_f[2*i],
                                2*i+1,sub_grad_f[2*i+1]);
                    }
                    printf("\n");
                    getchar();
                }

                part_UpdatePotentialGrid(x,x_bak,right_mId,right_mId.size());
                //                   UpdatePotentialGrid(x);    //time-consuming  0507

                for(int i = 0; i < right_mId.size();i++){
                    x_bak[ 2*right_mId[i] ]   = x[ 2*right_mId[i] ];
                    x_bak[ 2*right_mId[i]+1 ] = x[ 2*right_mId[i]+1 ];
                }
            }// end optimize one group net

            UpdateBlockPosition(x);
            UpdateDensityGrid(n,x);
            if(!param.bShow){
                maxDen = GetMaxDensity();
                totalOverDen = GetTotalOverDensity();
                totalOverDenLB = GetTotalOverDensityLB();
                totalOverPotential = GetTotalOverPotential();

                printf( " %d-%6d HPWL= %.0f\tDen= %.2f %.4f %.4f %.4f  LTime= %.1fm Dcost= %4.1f%% WireW= %.0f \n",
                        currentLevel, i_nets, m_pDB->CalcHPWL(),
                        maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                        double(seconds()-time_start)/60.0,
                        density*_weightDensity /objValue * 100.0,
                        _weightWire
                        );
            }

        }// end traversal nets
        sg_eval_f(x,obj2);
        if(obj2 < obj1)
            printf("down!down!down!down!down!down!down!down!down!down!down!down!down!down!down!down!down!down!down!down!\n");
        else
            printf("up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!up!\n");
        int merge_step = 1;
        if(ite % merge_step == 0){//
            double SL_step = 0;
            AdjustForce(n,x,SL_grad_f);
            LineSearch(n, x, SL_grad_f, SL_step);
            for(unsigned int i = 0; i < m_pDB->m_modules.size();i++){
                x[2*i]   += SL_step*SL_grad_f[2*i];
                x[2*i + 1] += SL_step*SL_grad_f[2*i + 1];
                lastPosition[i] = x[i];
            }

            BoundX(n, x, x_l, x_u);
            UpdateBlockPosition(x);
            UpdateDensityGrid(n,x);
            UpdatePotentialGrid(x);
        }




        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();
        if(param.bShow){
            printf( " %d-%2d HPWL= %.0f\tDen= %.2f %.4f %.4f %.4f  LTime= %.1fm Dcost= %4.1f%% WireW= %.0f \n",
                    currentLevel, ite+1, m_pDB->CalcHPWL(),
                    maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                    double(seconds()-time_start)/60.0,
                    density*_weightDensity /objValue * 100.0,
                    _weightWire
                    );
        }
#if 1
        if( param.bShow )
        {
            sprintf( filename, "SLfig%d-%d.plt", currentLevel, ite );
            m_pDB->OutputGnuplotFigure( filename, false, false );
        }

#endif
        if(ite % checkStep == 0 ){//record objValue per checkStep
            last_objValue = now_objValue;
            sg_eval_f(x, objValue);
            now_objValue = objValue;

            //                UpdateBlockPosition(x);// update modules positions per checkStep
            //                last_wirelength = now_wirelength;
            //                now_wirelength = m_pDB->CalcHPWL();

        }
        if( ite % checkStep == 0 )
        {
            printf( "." );
            fflush( stdout );

            if( ite % (2 * checkStep) == 0 )
            {
                UpdateBlockPosition( x );   // update to placeDB
                if( m_pDB->CalcHPWL() > bestLegalWL )   // gWL > LAL-WL
                {
                    printf( "BestLegalWL Break!\n" );
                    fflush( stdout );
                    break;
                }
            }

            UpdateDensityGrid( n, x );  // find the exact bin density
            totalOverDen = GetTotalOverDensity();
            totalOverDenLB = GetTotalOverDensityLB();
            totalOverPotential = GetTotalOverPotential();

            lastTotalOver = over;
            over = min( totalOverPotential, totalOverDen ); // TEST

            if( !startDecreasing
                    && over < lastTotalOver
                    && ite >= 6 )
            {
                printf( ">>" );
                fflush( stdout );
                startDecreasing = true;
            }

            // 2005-03-11: meet the constraint
            if( startDecreasing && over < target_density )
                break;

            // Cannot further improve
            if( now_objValue >= param.precision * last_objValue )
            {
                printf("Cannot further improve = _ =\n");
                for (int i = 0; i < n; i++) {
                    x[i] = lastPosition[i];
                }
                break;
            }

        } // end of ite == checkStep
        if(ite%1 == 0) {
            for (int i = 0; i < n; i++) {
                lastPosition[i] = x[i];
            }
        }

        //        }// end while true

        UpdateDensityGrid(n,x);
        UpdateBlockPosition(x);

        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();

        if (ite >= 1 && m_lookAheadLegalization && over < target_density + 0.1 )
        {
            UpdateBlockPosition(x);
            double hpwl = m_pDB->CalcHPWL();
            if (hpwl > bestLegalWL)
            {
                printf("Stop. Good enough.\n");
                break;
            }

            lookAheadLegalCount++;
            double oldWL = hpwl;
            CTetrisLegal legal(*m_pDB);

            double scale = 0.85;
            if (givenTargetUtil < 1.0 && givenTargetUtil > 0)
            {
                scale = 0.9;
            }

            double legalStart = seconds();
            bool   bLegal = legal.Solve( givenTargetUtil, false, false, scale);
            double legalTime = seconds() - legalStart;

            totalLegalTime += legalTime;

            if (param.bShow)
            {
                printf("LAL Time: %.2f\n", legalTime);
            }

            if (bLegal)
            {
                double WL = m_pDB->GetHPWLdensity(givenTargetUtil);
                if (param.bShow)
                {
                    m_pDB->ShowDensityInfo();
                }

                if (WL < bestLegalWL)
                {
                    LALnoGoodCount = 0;
                    if (param.bShow)
                    {
                        printf("SAVE BEST! (HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n",
                               m_pDB->GetHPWLp2p(), WL, (WL - oldWL) / oldWL * 100);
                    }
                    bestLegalWL = WL;
                    hasBestLegalSol = true;

                    for (int i = 0; i < (int) m_pDB->m_modules.size(); i++)
                    {
                        xBest[2 * i]     = m_pDB->m_modules[i].m_cx;
                        xBest[2 * i + 1] = m_pDB->m_modules[i].m_cy;
                        //                        x[2 * i]     = m_pDB->m_modules[i].m_cx;
                        //                        x[2 * i + 1] = m_pDB->m_modules[i].m_cy;

                    }
                } else {
                    if( param.bShow )
                        printf( "(HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n", m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    if( (WL-oldWL)/oldWL < 0.01 )
                    {
                        if( param.bShow )
                            printf( "Stop. Good enough\n" );
                        break;
                    }
                    LALnoGoodCount++;
                    if( LALnoGoodCount >= 2 )
                        break;
                }
            }
        } //end of m_lookAheadLegalization
        if (ite >= 1)
        {
            if ( m_lookAheadLegalization && startDecreasing && over < target_density + 0.01)
            {
                printf("Meet constraint!\n");
                break;
            }
            if (!m_lookAheadLegalization && startDecreasing && over < target_density + 0.03)
            {
                printf("Meet constraint!\n");
                break;
            }

            if (ite > 3 && totalOverPotential > lastTotalOverPotential && totalOverPotential < 1.4)
            {
                printf("Cannot further reduce!\n");
                break;
            }
        }
        // reduce _weightwire instead of _weightDensity  @hdgao
        /*
        if (over > 1.5)
        {
            _weightWire /= 1.5;
        } else if (over > 1.1) {
            _weightWire /= 1.5;
        } else {
            _weightWire /= 1.3;
        }
*/
        //        _weightDensity *=2;  //
        _weightWire /=2;

        lastTotalOverPotential = totalOverPotential;

    }// end of "for (ite = 0; ite < maxIte; ite++)"

    if (hasBestLegalSol)
    {
        memcpy(x, xBest, sizeof(double) * n);
    }

    UpdateBlockPosition(x);

    if( lookAheadLegalCount > 0 && param.bShow )
    {
        printf( "LAL: Total Count: %d\n", lookAheadLegalCount );
        printf( "LAL: Total CPU: %.2f\n", totalLegalTime );
        sprintf( filename, "util-global.dat" );
        CPlaceBin placeBin( *m_pDB );
        placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
        placeBin.OutputBinUtil( filename );
    }

    delete [] lastPosition;

    delete [] sub_grad_f;
    delete [] last_sub_grad_f;
    delete [] real_last_grad_f;
    delete [] old_grad_f;
    delete [] SL_grad_f;
    return hasBestLegalSol;
}
bool   MyNLP::SLSolve_smooth_bak(double wWire, double target_density, int currentLevel)
{
    double givenTargetUtil = m_targetUtil;    // m_targetUtil = -1
    m_currentStep = param.step;

    m_targetUtil += 0.05;
    if (m_targetUtil > 1.0)
        m_targetUtil = 1.0;

    double time_start = seconds();
    char filename[100];    // for gnuplot

    int n = 2 * m_pDB->m_modules.size();  //for the total num of modules

    // calculate the goal of the utilization
    double designUtil = m_pDB->m_totalModuleArea / m_pDB->m_totalFreeSpace;

    if (param.bShow)
        printf( "hdgao INFO: Design utilization: %f\n", designUtil);

    if (m_targetUtil > 0)
    {  // has give a utilization
        double lowest = designUtil + 0.05;
        if (m_targetUtil < lowest)
        {
            if (param.bShow)
            {
                printf("hdgao WARNING: Target utilization (%f) is too low\n", m_targetUtil);
                printf("hdgao          Set target utilization to %f \n", lowest);
            }
            m_targetUtil = lowest;
        }
    } else { // no given utilization
        printf("flpeng WARNING: No given target utilization. Distribute blocks evenly.\n");
        m_targetUtil = designUtil + 0.05;
        if (m_targetUtil > 1.0)
            m_targetUtil = 1.0;
    }

    if (param.bShow)
        printf("flpeng DBIN: Target utilization: %f\n", m_targetUtil);



    //flpengchange
    double* lastPosition = new double[n];
    memset(lastPosition, 0, sizeof(double) * n);


    //calculate the subgrad of f and the direction and allocate memeory
    double* sub_grad_f = new double [ n];
    double* last_sub_grad_f = new double [ n ]; // for computing CG-direction;
    double* real_last_grad_f = new double [ n];
    double* old_grad_f = new double[n];
    memset( sub_grad_f, 0, sizeof(double)*n );
    memset( last_sub_grad_f, 0, sizeof(double)*n );
    memset(real_last_grad_f, 0, sizeof(double) * n);
    memset(old_grad_f, 0, sizeof(double) *n);

    double* SL_grad_f = new double [n];
    memset( SL_grad_f, 0, sizeof(double)*n );

    CreatePotentialGrid();   // create potential grid according to "m_potentialGridSize"
    int densityGridSize = 10;	// 1% chip area
    //int densityGridSize = m_potentialGridSize / 3;	    // not good in big3
    CreateDensityGrid( densityGridSize );	// real density: use 1% area
    //flpeng add so far the x has been calculated by qsolve or the last uncoarsening iteration
    UpdateDensityGridSpace( n, x );
    UpdatePotentialGridBase( x );		// init exp potential for each bin, also update ExpBin


#if 1
    // gaussian smoothing for base potential
    GaussianSmooth smooth;
    int r = m_smoothR;
    smooth.Gaussian2D( r, 6*r+1 );
    smooth.Smooth( m_basePotential );
    m_basePotentialOri = m_basePotential;
    sprintf( filename, "gbase%d.dat", currentLevel );
    OutputPotentialGrid( filename );

    cout << "\nDebug: m_smoothDelta is : " << m_smoothDelta << endl;
    //  cout << "\nDebug: m_smoothDelta is : " << m_smoothDelta << endl;

    // TEST
    if( m_smoothDelta == 1 )
    {
        if( param.bShow )
        {
            sprintf( filename, "gbase%d-more.dat", currentLevel );
            printf( "generate %s...\n", filename );
            fflush( stdout );
        }

        vector< vector< double > > moreSmooth = m_basePotential;
        r = m_smoothR * 6;
        int kernel_size = 5*r;
        if( kernel_size % 2 == 0 )
            kernel_size++;
        smooth.Gaussian2D( r, kernel_size );
        smooth.Smooth( moreSmooth );

        if( param.bShow )
        {
            swap( moreSmooth, m_basePotential );
            OutputPotentialGrid( filename );
            swap( moreSmooth, m_basePotential );
        }

        // merge base and moreSmooth
        double binArea = m_potentialGridWidth * m_potentialGridHeight;
        double halfBinArea = binArea / 2;
        int changeCount = 0;
        for( unsigned int i=0; i<moreSmooth.size(); i++ )
        {
            for( unsigned int j=0; j<moreSmooth[i].size(); j++ )
            {
                double free = binArea - m_basePotential[i][j];
                if( free < 1e-4 )	// no space
                {
                    if( moreSmooth[i][j] > halfBinArea )
                    {
                        m_basePotential[i][j] += moreSmooth[i][j] - halfBinArea;
                        changeCount++;
                    }
                }
            }
        }

        if( param.bShow )
        {
            printf( "change %d\n", changeCount );
            sprintf( filename, "gbase%d-more-merge.dat", currentLevel );
            OutputPotentialGrid( filename );
        }
    }


    if( m_smoothDelta > 1.0 )
        SmoothPotentialBase( double(m_smoothDelta) );   // also update ExpBin

    UpdateExpBinPotential( m_targetUtil );

    if( param.bShow )
    {
        sprintf( filename, "base%d.dat", currentLevel );
        OutputPotentialGrid( filename );
    }
#endif
    double objValue;

    assert( m_targetUtil > 0 );

    cout << "\nDebug: the inAlpha: " << _alpha << '\t' << "  the scale: " << m_posScale << endl;
    // wirelength
    UpdateExpValueForEachCell( n, x, _expX, _alpha );
    UpdateExpValueForEachPin( n, x, _expPins, _alpha );
    UpdateNetsSumExp( x, _expX );
    totalWL = GetWL( n, x, _expX, _alpha );

    // density
    UpdatePotentialGrid( x );
    UpdateDensityGrid( n, x );
    density = GetDensityPanelty();


    // 2006-02-22 weight (APlace ICCAD05)
    _weightWire = 1.0;
    eval_grad_f( n, x, _expX, true, sub_grad_f );
    double totalWireGradient = 0;
    double totalPotentialGradient = 0;

    // TODO: truncation?
    AdjustForce( n, x, grad_wire, grad_potential );

    for( int i=0; i<n; i++ )
    {
        totalWireGradient      += fabs( grad_wire[i] );
        totalPotentialGradient += fabs( grad_potential[i] );
    }


    _weightWire = 1.0;

    density = GetDensityPanelty();// to evaluate overlap area  @2704hdgao
    sg_eval_f(x, objValue);
    sg_eval_grad_f(n, x, sub_grad_f);
    AdjustForce(n, x, sub_grad_wire, sub_grad_potential); // ? modify sub_grad_wire and sub_grad_potential

    double totalSubWireGradient = 0.0;
    double totalSubPotentialGradient = 0.0;
    for (int i = 0; i < n; i++) {
        totalSubWireGradient      += fabs(sub_grad_wire[i]);
        totalSubPotentialGradient += fabs(sub_grad_potential[i]);
    }

    _weightDensity = 1.0;
    _weightWire = wWire * totalSubPotentialGradient / totalSubWireGradient ;//* param.sg_namdaIncFactor ;

    int      maxIte = 50;   // max iterator
    double   beta;   // for the CG method
    double nnbReal = GetNonZeroDensityGridPercent();
    double maxDen = GetMaxDensity();
    double totalOverDen = GetTotalOverDensity();
    double totalOverDenLB = GetTotalOverDensityLB();
    double totalOverPotential = GetTotalOverPotential();

    if (param.bShow)
    {
        printf(" %d-%2d HPWL= %.0f\tDen= %.2f %.2f %.2f %.2f NNB= %.2f Dcost= %4.1f%%  WireW= %.0f",
               currentLevel, m_ite, m_pDB->CalcHPWL(), maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
               nnbReal, density * _weightDensity / objValue * 100.0, _weightWire
               );
    } else {
        printf( " %d-%2d HPWL= %.0f \t", currentLevel, m_ite, m_pDB->CalcHPWL() );
    }
    fflush (stdout);

    if( param.bShow )
    {
        sprintf( filename, "SLfig%d-%d.plt", currentLevel, m_ite );
        m_pDB->OutputGnuplotFigure( filename, false, false );
    }


    double lastTotalOver = 0.0;
    double lastTotalOverPotential = DBL_MAX;
    double over = totalOverDen;
    //    int    totalIte = 0;
    bool   hasBestLegalSol = false;
    double bestLegalWL = DBL_MAX;
    int    lookAheadLegalCount = 0;
    double totalLegalTime = 0.0;
    bool   startDecreasing = false;
    int checkStep = 1;
    int LALnoGoodCount = 0;

    srand(time(NULL));
    double ssModify = 1.0;
    stepSize = 1; // 1 - > 2





    printf("\n======================================================================================\n");
    printf("cg begin\n");
    for(int ite = 0;ite < maxIte;ite++){
        //        while(true){
        //            int update_counts = 0;
        m_currentStep = param.sg_stepsize;
        memset( SL_grad_f, 0, sizeof(double)*n );

        double last_objValue = DBL_MAX;
        double now_objValue = DBL_MAX;

        //            double last_wirelength = DBL_MAX;
        //           double now_wirelength  = DBL_MAX;
        //            double best_wirelenth = DBL_MAX;

        //  update modules positions
        int num_nets = 60;
        bool nets_flag = true;
        unsigned int total_nets = m_pDB->m_nets.size();
        for(unsigned int i_nets = 0; nets_flag ;i_nets += num_nets){

            if(i_nets + num_nets >= total_nets){
                num_nets = total_nets - i_nets + 1;
                nets_flag = false;
            }

            int current_module_size = 0;
            for(int i = 0;i<num_nets;i++){
                current_module_size += m_pDB->m_nets[ i+i_nets ].size();
            }

            memset( sub_grad_f, 0, sizeof(double)*2*current_module_size);
            memset( last_sub_grad_f, 0, sizeof(double)*2*current_module_size);
            memset(real_last_grad_f, 0, sizeof(double) * 2 * current_module_size);
            memset(old_grad_f, 0, sizeof(double) * 2 * current_module_size);

            //                part_eval_f(x, objValue,i_nets);// get objValue
            smooth_part_eval_grad_f(x,sub_grad_f,i_nets,num_nets); // get grad
            AdjustForce(2*current_module_size,x,sub_grad_f);    // m_gridDensity is  irrelevant to sub_grad_f

            for(int i_cg = 0;i_cg < 5;i_cg++){// 5 -> 3
                swap(real_last_grad_f, old_grad_f);
                swap(last_sub_grad_f, sub_grad_f);
                smooth_part_eval_grad_f(x, sub_grad_f,i_nets,num_nets);
                AdjustForce(2*current_module_size,x,sub_grad_f);
                //calculate d_k
                if (i_cg == 0){// first
                    for (int i = 0; i < 2*current_module_size; i++){
                        old_grad_f[i] = sub_grad_f[i];
                        sub_grad_f[i] = -sub_grad_f[i];
                    }
                }
                else {
                    // conjugate sub gradient direction
                    sg_findBeta(2*current_module_size, sub_grad_f, last_sub_grad_f, real_last_grad_f, beta);
                    // FindBeta(n, sub_grad_f, last_sub_grad_f, beta);
                    for (int i = 0; i < 2*current_module_size; i++)
                    {
                        old_grad_f[i] = sub_grad_f[i];
                        sub_grad_f[i] = -sub_grad_f[i] + beta * last_sub_grad_f[i];
                    }
                }

                //calculate a_k , the step size

                LineSearch(2*current_module_size, x, sub_grad_f, stepSize);

                for(int j = 0,k = 0;j < num_nets;j++){
                    for(int i = 0;i < m_pDB->m_nets[i_nets + j].size();i++){
                        int pinId =  m_pDB->m_nets[i_nets + j][i];
                        int moduleId =  m_pDB->m_pins[pinId].moduleId;
                        x[2*moduleId]   += sub_grad_f[2*k]   * stepSize;
                        x[2*moduleId+1] += sub_grad_f[2*k+1] * stepSize;
                        // record the final part_grad
                        SL_grad_f[2*moduleId]     += sub_grad_f[2*k];
                        SL_grad_f[2*moduleId+1]   += sub_grad_f[2*k+1];
                        k++;
                    }
                }
                BoundX(n, x, x_l, x_u);
                double time_used = seconds();
                UpdateExpValueForEachCell( n, x, _expX, _alpha );
                UpdateExpValueForEachPin( n, x, _expPins, _alpha );
                UpdateNetsSumExp( x, _expX );
                time_grad_wl += seconds() - time_used;
                UpdatePotentialGrid( x );

            }// end optimize one group net

            UpdateBlockPosition(x);
            UpdateDensityGrid(n,x);
            if(!param.bShow){
                maxDen = GetMaxDensity();
                totalOverDen = GetTotalOverDensity();
                totalOverDenLB = GetTotalOverDensityLB();
                totalOverPotential = GetTotalOverPotential();

                printf( " %d-%6d HPWL= %.0f\tDen= %.2f %.4f %.4f %.4f  LTime= %.1fm Dcost= %4.1f%% WireW= %.0f \n",
                        currentLevel, i_nets, m_pDB->CalcHPWL(),
                        maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                        double(seconds()-time_start)/60.0,
                        density*_weightDensity /objValue * 100.0,
                        _weightWire
                        );
            }

        }// end traversal nets
#if 0
        int merge_step = 1;
        if(ite % merge_step == 0){//
            double SL_step = 0;
            AdjustForce(n,x,SL_grad_f);
            LineSearch(n, x, SL_grad_f, SL_step);
            for(unsigned int i = 0; i < m_pDB->m_modules.size();i++){
                x[2*i]   += SL_step*SL_grad_f[2*i];
                x[2*i + 1] += SL_step*SL_grad_f[2*i + 1];
                lastPosition[2*i] = x[2*i];
                lastPosition[2*i+1] = x[2*i+1];
            }
            BoundX( n, x, x_l, x_u );
            UpdateExpValueForEachCell( n, x, _expX, _alpha );
            UpdateExpValueForEachPin( n, x, _expPins, _alpha );
            UpdateNetsSumExp( x, _expX );
            UpdateBlockPosition(x);
            UpdateDensityGrid(n,x);
            UpdatePotentialGrid(x);

        }
#endif
        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();
        if(param.bShow){
            printf( " %d-%2d HPWL= %.0f\tDen= %.2f %.4f %.4f %.4f  LTime= %.1fm Dcost= %4.1f%% WireW= %.0f \n",
                    currentLevel, ite+1, m_pDB->CalcHPWL(),
                    maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                    double(seconds()-time_start)/60.0,
                    density*_weightDensity /objValue * 100.0,
                    _weightWire
                    );
        }
#if 1
        if( param.bShow )
        {
            sprintf( filename, "SLfig%d-%d.plt", currentLevel, ite );
            m_pDB->OutputGnuplotFigure( filename, false, false );
        }

#endif
        if(ite % checkStep == 0 ){//record objValue per checkStep
            last_objValue = now_objValue;
            eval_f( n, x, _expX, true, objValue );
            now_objValue = objValue;

            //                UpdateBlockPosition(x);// update modules positions per checkStep
            //                last_wirelength = now_wirelength;
            //                now_wirelength = m_pDB->CalcHPWL();

        }
        if( ite % checkStep == 0 )
        {
            printf( "." );
            fflush( stdout );

            if( ite % (2 * checkStep) == 0 )
            {
                UpdateBlockPosition( x );   // update to placeDB
                if( m_pDB->CalcHPWL() > bestLegalWL )   // gWL > LAL-WL
                {
                    printf( "BestLegalWL Break!\n" );
                    fflush( stdout );
                    break;
                }
            }

            UpdateDensityGrid( n, x );  // find the exact bin density
            totalOverDen = GetTotalOverDensity();
            totalOverDenLB = GetTotalOverDensityLB();
            totalOverPotential = GetTotalOverPotential();

            lastTotalOver = over;
            over = min( totalOverPotential, totalOverDen ); // TEST

            if( !startDecreasing
                    && over < lastTotalOver
                    && ite >= 6 )
            {
                printf( ">>" );
                fflush( stdout );
                startDecreasing = true;
            }

            // 2005-03-11: meet the constraint
            if( startDecreasing && over < target_density )
                break;

            // Cannot further improve
            if( now_objValue >= param.precision * last_objValue )
            {
                printf("Cannot further improve = _ =\n");
                for (int i = 0; i < n; i++) {
                    x[i] = lastPosition[i];
                }
                break;
            }

        } // end of ite == checkStep
        if(ite%1 == 0) {
            for (int i = 0; i < n; i++) {
                lastPosition[i] = x[i];
            }
        }

        //        }// end while true

        UpdateDensityGrid(n,x);
        UpdateBlockPosition(x);

        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();

        if (ite >= 1 && m_lookAheadLegalization && over < target_density + 0.1 )
        {
            UpdateBlockPosition(x);
            double hpwl = m_pDB->CalcHPWL();
            if (hpwl > bestLegalWL)
            {
                printf("Stop. Good enough.\n");
                break;
            }

            lookAheadLegalCount++;
            double oldWL = hpwl;
            CTetrisLegal legal(*m_pDB);

            double scale = 0.85;
            if (givenTargetUtil < 1.0 && givenTargetUtil > 0)
            {
                scale = 0.9;
            }

            double legalStart = seconds();
            bool   bLegal = legal.Solve( givenTargetUtil, false, false, scale);
            double legalTime = seconds() - legalStart;

            totalLegalTime += legalTime;

            if (param.bShow)
            {
                printf("LAL Time: %.2f\n", legalTime);
            }

            if (bLegal)
            {
                double WL = m_pDB->GetHPWLdensity(givenTargetUtil);
                if (param.bShow)
                {
                    m_pDB->ShowDensityInfo();
                }

                if (WL < bestLegalWL)
                {
                    LALnoGoodCount = 0;
                    if (param.bShow)
                    {
                        printf("SAVE BEST! (HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n",
                               m_pDB->GetHPWLp2p(), WL, (WL - oldWL) / oldWL * 100);
                    }
                    bestLegalWL = WL;
                    hasBestLegalSol = true;

                    for (int i = 0; i < (int) m_pDB->m_modules.size(); i++)
                    {
                        xBest[2 * i]     = m_pDB->m_modules[i].m_cx;
                        xBest[2 * i + 1] = m_pDB->m_modules[i].m_cy;
                        //                        x[2 * i]     = m_pDB->m_modules[i].m_cx;
                        //                        x[2 * i + 1] = m_pDB->m_modules[i].m_cy;

                    }
                } else {
                    if( param.bShow )
                        printf( "(HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n", m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    if( (WL-oldWL)/oldWL < 0.01 )
                    {
                        if( param.bShow )
                            printf( "Stop. Good enough\n" );
                        break;
                    }
                    LALnoGoodCount++;
                    if( LALnoGoodCount >= 2 )
                        break;
                }
            }
        } //end of m_lookAheadLegalization
        if (ite >= 1)
        {
            if ( m_lookAheadLegalization && startDecreasing && over < target_density + 0.01)
            {
                printf("Meet constraint!\n");
                break;
            }
            if (!m_lookAheadLegalization && startDecreasing && over < target_density + 0.03)
            {
                printf("Meet constraint!\n");
                break;
            }

            if (ite > 3 && totalOverPotential > lastTotalOverPotential && totalOverPotential < 1.4)
            {
                printf("Cannot further reduce!\n");
                break;
            }
        }
        // reduce _weightwire instead of _weightDensity  @hdgao
        /*
        if (over > 1.5)
        {
            _weightWire /= 1.5;
        } else if (over > 1.1) {
            _weightWire /= 1.5;
        } else {
            _weightWire /= 1.3;
        }
*/
        _weightWire /= m_incFactor;   //  _weightDensity *=2;
        //       _weightDensity *=2;
        lastTotalOverPotential = totalOverPotential;

    }// end of "for (ite = 0; ite < maxIte; ite++)"

    if (hasBestLegalSol)
    {
        memcpy(x, xBest, sizeof(double) * n);
    }

    UpdateBlockPosition(x);

    if( lookAheadLegalCount > 0 && param.bShow )
    {
        printf( "LAL: Total Count: %d\n", lookAheadLegalCount );
        printf( "LAL: Total CPU: %.2f\n", totalLegalTime );
        sprintf( filename, "util-global.dat" );
        CPlaceBin placeBin( *m_pDB );
        placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
        placeBin.OutputBinUtil( filename );
    }

    delete [] lastPosition;

    delete [] sub_grad_f;
    delete [] last_sub_grad_f;
    delete [] real_last_grad_f;
    delete [] old_grad_f;
    delete [] SL_grad_f;
    return hasBestLegalSol;
}

bool   MyNLP::SLSolve_smooth(double wWire, double target_density, int currentLevel)
{
    double givenTargetUtil = m_targetUtil;    // m_targetUtil = -1
    m_currentStep = param.step;

    m_targetUtil += 0.05;
    if (m_targetUtil > 1.0)
        m_targetUtil = 1.0;

    double time_start = seconds();
    char filename[100];    // for gnuplot

    int n = 2 * m_pDB->m_modules.size();  //for the total num of modules

    // calculate the goal of the utilization
    double designUtil = m_pDB->m_totalModuleArea / m_pDB->m_totalFreeSpace;

    if (param.bShow)
        printf( "hdgao INFO: Design utilization: %f\n", designUtil);

    if (m_targetUtil > 0)
    {  // has give a utilization
        double lowest = designUtil + 0.05;
        if (m_targetUtil < lowest)
        {
            if (param.bShow)
            {
                printf("hdgao WARNING: Target utilization (%f) is too low\n", m_targetUtil);
                printf("hdgao          Set target utilization to %f \n", lowest);
            }
            m_targetUtil = lowest;
        }
    } else { // no given utilization
        printf("hdgao WARNING: No given target utilization. Distribute blocks evenly.\n");
        m_targetUtil = designUtil + 0.05;
        if (m_targetUtil > 1.0)
            m_targetUtil = 1.0;
    }

    if (param.bShow)
        printf("hdgao DBIN: Target utilization: %f\n", m_targetUtil);



    //flpengchange
    double* lastPosition = new double[n];
    memset(lastPosition, 0, sizeof(double) * n);


    //calculate the subgrad of f and the direction and allocate memeory
    double* p_grad_f      = new double[ 2*n ];
    double* last_p_grad_f = new double[ 2*n ];
    memset(p_grad_f,     0,sizeof(double)*n);
    memset(last_p_grad_f,0,sizeof(double)*n);

    double* SL_grad_f = new double [n];
    memset( SL_grad_f, 0, sizeof(double)*n );

    CreatePotentialGrid();   // create potential grid according to "m_potentialGridSize"
    int densityGridSize = 10;	// 1% chip area
    //int densityGridSize = m_potentialGridSize / 3;	    // not good in big3
    CreateDensityGrid( densityGridSize );	// real density: use 1% area
    //flpeng add so far the x has been calculated by qsolve or the last uncoarsening iteration
    UpdateDensityGridSpace( n, x );
    UpdatePotentialGridBase( x );		// init exp potential for each bin, also update ExpBin


#if 1
    // gaussian smoothing for base potential
    GaussianSmooth smooth;
    int r = m_smoothR;
    smooth.Gaussian2D( r, 6*r+1 );
    smooth.Smooth( m_basePotential );
    m_basePotentialOri = m_basePotential;
    sprintf( filename, "gbase%d.dat", currentLevel );
    OutputPotentialGrid( filename );
#endif

    // TEST
    if( m_smoothDelta == 1 )
    {
        if( param.bShow )
        {
            sprintf( filename, "gbase%d-more.dat", currentLevel );
            printf( "generate %s...\n", filename );
            fflush( stdout );
        }

        vector< vector< double > > moreSmooth = m_basePotential;
        r = m_smoothR * 6;
        int kernel_size = 5*r;
        if( kernel_size % 2 == 0 )
            kernel_size++;
        smooth.Gaussian2D( r, kernel_size );
        smooth.Smooth( moreSmooth );

        if( param.bShow )
        {
            swap( moreSmooth, m_basePotential );
            OutputPotentialGrid( filename );
            swap( moreSmooth, m_basePotential );
        }

        // merge base and moreSmooth
        double binArea = m_potentialGridWidth * m_potentialGridHeight;
        double halfBinArea = binArea / 2;
        int changeCount = 0;
        for( unsigned int i=0; i<moreSmooth.size(); i++ )
        {
            for( unsigned int j=0; j<moreSmooth[i].size(); j++ )
            {
                double free = binArea - m_basePotential[i][j];
                if( free < 1e-4 )	// no space
                {
                    if( moreSmooth[i][j] > halfBinArea )
                    {
                        m_basePotential[i][j] += moreSmooth[i][j] - halfBinArea;
                        changeCount++;
                    }
                }
            }
        }

        if( param.bShow )
        {
            printf( "change %d\n", changeCount );
            sprintf( filename, "gbase%d-more-merge.dat", currentLevel );
            OutputPotentialGrid( filename );
        }
    }


    if( m_smoothDelta > 1.0 )
        SmoothPotentialBase( double(m_smoothDelta) );   // also update ExpBin

    UpdateExpBinPotential( m_targetUtil );

    if( param.bShow )
    {
        sprintf( filename, "base%d.dat", currentLevel );
        OutputPotentialGrid( filename );
    }

    double objValue;

    assert( m_targetUtil > 0 );

    // wirelength
    UpdateExpValueForEachCell( n, x, _expX, _alpha );
    UpdateExpValueForEachPin( n, x, _expPins, _alpha );
    UpdateNetsSumExp( x, _expX );
    totalWL = GetWL( n, x, _expX, _alpha );

    // density
    UpdatePotentialGrid( x );
    UpdateDensityGrid( n, x );
    density = GetDensityPanelty();


    // 2006-02-22 weight (APlace ICCAD05)
    _weightWire = 1.0;
    eval_grad_f( n, x, _expX, true, p_grad_f );
    double totalWireGradient = 0;
    double totalPotentialGradient = 0;

    // TODO: truncation?
    AdjustForce( n, x, grad_wire, grad_potential );

    for( int i=0; i<n; i++ )
    {
        totalWireGradient      += fabs( grad_wire[i] );
        totalPotentialGradient += fabs( grad_potential[i] );
    }


    _weightDensity = 1.0*1.02;
    _weightWire = wWire *  totalPotentialGradient/totalWireGradient;//*param.sg_namdaIncFactor;//* param.sg_namdaIncFactor ;

    int      maxIte = 50;   // max iterator
    //    double   beta;   // for the CG method
    double obj_value;
    eval_f(n,x,_expX,true,obj_value);
    eval_grad_f(n,x,_expX,true,p_grad_f);

    double nnbReal = GetNonZeroDensityGridPercent();
    UpdateDensityGrid( n, x );
    double maxDen = GetMaxDensity();
    double totalOverDen = GetTotalOverDensity();
    double totalOverDenLB = GetTotalOverDensityLB();
    double totalOverPotential = GetTotalOverPotential();

    if (param.bShow)
    {
        printf(" %d-%2d HPWL= %.0f\tDen= %.2f %.2f %.2f %.2f NNB= %.2f Dcost= %4.1f%%  WireW= %.0f",
               currentLevel, m_ite, m_pDB->CalcHPWL(), maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
               nnbReal, density * _weightDensity / obj_value * 100.0, _weightWire
               );
    } else {
        printf( " %d-%2d HPWL= %.0f \t", currentLevel, m_ite, m_pDB->CalcHPWL() );
    }
    fflush (stdout);

    if( param.bShow )
    {
        sprintf( filename, "SLfig%d-%d.plt", currentLevel, m_ite );
        m_pDB->OutputGnuplotFigure( filename, false, false );
    }


    double lastTotalOver = 0.0;
    double lastTotalOverPotential = DBL_MAX;
    double over = totalOverDen;
    //    int    totalIte = 0;
    bool   hasBestLegalSol = false;
    double bestLegalWL = DBL_MAX;
    int    lookAheadLegalCount = 0;
    double totalLegalTime = 0.0;

    bool   startDecreasing = false;

    int checkStep = 5;

    int LALnoGoodCount = 0;
    int totalIte = 0;

    int bucket_size = 10;
    vector< vector< vector<int> > > bucket;
    test();
    if(1 && currentLevel != 1){
#if 0
        NCA_deeper();
        nets_cluster_modules_bounds.resize( nets_cluster.size() );
        for(int i = 0;i < nets_cluster.size();i++){
            getBlackModule(i);
        }
#endif

        double bucket_width  = (m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left)/bucket_size;
        double bucket_height = (m_pDB->m_coreRgn.top - m_pDB->m_coreRgn.bottom)/bucket_size;



        bucket.resize(bucket_size);
        for(int i = 0;i < bucket_size;i++)
            bucket[i].resize(bucket_size);

        //        for(int i = 0;i < bucket_size;i++)
        //            for(int j = 0;j < bucket_size;j++)
        //                bucket[i][j].clear();

        for(unsigned int i = 0;i < m_pDB->m_nets.size();i++){
            double c_x,c_y;
            get_center_position(x,c_x,c_y,i,0);//0: ortho_center,1: mean_center
            int n_x = (c_x - m_pDB->m_coreRgn.left)/(bucket_width );
            int n_y = (c_y - m_pDB->m_coreRgn.bottom)/(bucket_height);
            //            if(n_x == bucket_size)
            //                n_x--;
            //            if(n_y == bucket_size)
            //                n_y--;
            bucket[n_x][n_y].push_back((int)i);
        }


        if(param.bShow){
            printf("\n");
            int t = 0;
            for(int i = 0;i < bucket_size;i++){
                for(int j = 0;j < bucket_size;j++){
                    printf("%4d\t",bucket[i][j].size());
                    t += bucket[i][j].size();
                }
                printf("\n");
            }
            printf("net size  %d\n",t);
        }
        //        getchar();
    }
    else{//currentLevel = 1
        //whole net

        nets_cluster.resize(1);
        for(int i = 0;i < m_pDB->m_nets.size();i++)
            nets_cluster[0].push_back(i);
        nets_cluster_modules_bounds.resize( nets_cluster.size() );
        getBlackModule(0);
    }


    //   int cycle_limit = (currentLevel == 1) ? 100000:param.run_cycle;
    double* x_bak = new double[n];

    printf("\n======================================================================================\n");
    printf("cg begin\n");
    bool is_2D = false;
    for(int ite = 0;ite < maxIte;ite++){
        int innerIte = 0;
        double old_obj = DBL_MAX;
        double last_obj_value = DBL_MAX;


        m_currentStep = param.step;//param.sg_stepsize
        if(1 && currentLevel != 1)
            m_currentStep = 0.05;//param.run_cycle;
        for(int i = 0; i < n;i++){
            x_bak[i] = x[i];
        }

        while(true){
            innerIte++;
            //          if(innerIte > 1)  //test
            //             break;
            if( innerIte % checkStep == 0){
                old_obj = last_obj_value;
                eval_f(n,x,_expX,true,obj_value);
                last_obj_value = obj_value;
            }
            if( innerIte % checkStep == 0){
                printf(".");
                if( innerIte % (2*checkStep) == 0){
                    UpdateBlockPosition( x );
                    if( m_pDB->CalcHPWL() > bestLegalWL ){
                        printf("best legal WL\n");
                        break;
                    }
                }

                UpdateDensityGrid(n,x);
                totalOverDen = GetTotalOverDensity();
                totalOverDenLB = GetTotalOverDensityLB();
                totalOverPotential = GetTotalOverPotential();

                lastTotalOver = over;
                over = min( totalOverPotential, totalOverDen ); // TEST

                if( !startDecreasing
                        && over < lastTotalOver
                        && ite >= 1
                        && innerIte >= 6)
                {
                    printf( ">>" );
                    startDecreasing = true;
                }

                // 2005-03-11: meet the constraint
                if( startDecreasing && over < target_density )
                    break;
                if( obj_value >= param.precision * old_obj){
                    break;
                }
            }
            /*  pick random net */
            if(1 && currentLevel != 1){

                double pick_rate = 0.01;
                vector< vector<int> > temp_nets_cluster;
                nets_cluster.clear();
                nets_black_modules.clear();
                nets_cluster_modules_bounds.clear();

                temp_nets_cluster.resize(bucket_size * bucket_size);//temp_nets_cluster[0][0:99]

                for(int i = 0;i < bucket_size;i++)
                    for(int j = 0;j < bucket_size;j++){
                        int num_pick_net = pick_rate*bucket[i][j].size();
                        if(bucket[i][j].size() == 0) continue;
                        if(num_pick_net == 0) num_pick_net = 1;
                        vector<int> rand_vec = randVector(bucket[i][j].size());
                        for(int k = 0;k < num_pick_net;k++)
                            temp_nets_cluster[i*bucket_size + j].push_back(bucket[i][j][ rand_vec[k] ]);
                    }
                vector<int> rand_vec = randVector(bucket_size*bucket_size);
                for(int k = 0;k < bucket_size;k++){// this put all(100) regions into the nets_cluster
                    vector<int> temp;
                    for(int i = 0;i < bucket_size*1;i++)
                        temp.insert(temp.end(),temp_nets_cluster[ rand_vec[ k*bucket_size + i] ].begin(),
                                temp_nets_cluster[ rand_vec[ k*bucket_size + i] ].end() );
                    nets_cluster.push_back(temp);// nets_cluster[0][...]
                }



                nets_black_modules.resize(nets_cluster.size());// update nets_black_modules,nets_cluster_bounds according to nets_cluster
                nets_cluster_modules_bounds.resize(nets_cluster.size());
                for(unsigned int i = 0; i < nets_cluster.size();i++)
                    getBlackModule(i);
            }

            //  update modules positions
            /* if(currentLevel == 1)
                CG(this,n,x_bak,is_2D);
            else{// parallel
                do{
                    for(int i = 0;i < bucket_size;i++){
                        vector<int> specified_nets_cluster_Id;
                        specified_nets_cluster_Id.push_back(i);
                        CG(this,n,x_bak,is_2D,specified_nets_cluster_Id);
                    }
                }while(done == bucket_size);
            }*/


            unsigned int size_nets_cluster = nets_cluster.size();
            for(unsigned int i_nets_cluster = 0; i_nets_cluster < size_nets_cluster ;i_nets_cluster++){

                vector< double > last_grad_f,grad_f;
                last_grad_f.resize(n,0);
                grad_f.resize(n,0);
                vector<int> RegionModule;
                RegionModule.clear();
                //              RegionModule.assign(nets_black_modules[i_nets_cluster].begin(),nets_black_modules[i_nets_cluster].end());
                int cg_numIte = 5;
                for(int i_cg = 0;i_cg < cg_numIte;i_cg++){// 5 -> 3
                    last_grad_f.assign( grad_f.begin(),grad_f.end() );
                    if(i_cg == 0)  updateBlackBounds(i_nets_cluster);
                    smooth_part_eval_grad_f(x,nets_cluster_modules_bounds[i_nets_cluster],nets_black_modules[i_nets_cluster],
                                            RegionModule,grad_f); // need change

                    //                    double grad_square_sum = 0;
                    //                    for(int i = 0; i < RegionModule.size();i++){
                    //                        int moduleId = RegionModule[i];
                    //                        grad_square_sum += grad_f[2*moduleId] * grad_f[2*moduleId] + grad_f[2*moduleId+1] * grad_f[2*moduleId+1];
                    //                    }
                    //                    printf("%f\t",grad_square_sum);
                    //                    if(i_cg == cg_numIte - 1) printf("\n");


                    if(is_2D)
                        PartAdjustForce_2D(grad_f,RegionModule);
                    else
                        PartAdjustForce(grad_f,RegionModule);
                    //calculate d_k
                    double beta1,beta2,beta;
                    if(i_cg == 0){// first
                        for (unsigned int i = 0; i < RegionModule.size(); i++){
                            int moduleId = RegionModule[i];
                            grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ];
                            grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1];
                        }
                    }
                    else{
                        if(is_2D)
                            PartFindBeta_2D(grad_f,last_grad_f,RegionModule,beta1,beta2);
                        else{
                            PartFindBeta(grad_f,last_grad_f,RegionModule,beta);
                            beta1 = beta2 = beta;
                        }
                        for (unsigned int i = 0; i < RegionModule.size(); i++){
                            int moduleId = RegionModule[i];
                            grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ]    + beta1 * last_grad_f[ 2*moduleId ];
                            grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1] + beta2 * last_grad_f[ 2*moduleId + 1];
                        }
                    }

                    //calculate a_k , the step size
                    double stepSize1,stepSize2;
                    if(is_2D)
                        PartLineSearch_2D(grad_f,RegionModule,stepSize1,stepSize2);
                    else{
                        PartLineSearch(grad_f,RegionModule,stepSize);
                        stepSize1 = stepSize2 = stepSize;
                    }
                    for (int i = 0; i < RegionModule.size(); i++){
                        int moduleId = RegionModule[i];
                        x[2*moduleId]   += grad_f[2*moduleId]  * stepSize1;
                        x[2*moduleId+1] += grad_f[2*moduleId+1]* stepSize2;
                    }

                    BoundX(n, x, x_l, x_u);
                    double time_used = seconds();
                    vector<int> right_mId;
                    for(int i = 0; i < n/2;i++){
                        if( x_bak[2*i] != x[2*i] || x_bak[2*i + 1] != x[2*i + 1])
                            right_mId.push_back(i);
                    }
                    part_UpdatePotentialGrid(x,x_bak,right_mId,right_mId.size());
                    for(int i = 0; i < right_mId.size();i++){
                        x_bak[ 2*right_mId[i] ]   = x[ 2*right_mId[i] ];
                        x_bak[ 2*right_mId[i]+1 ] = x[ 2*right_mId[i]+1 ];
                    }

                    // time-consuming here

                    UpdateExpValueForEachCell( n, x, _expX, _alpha );
                    UpdateExpValueForEachPin( n, x, _expPins, _alpha );
                    UpdateNetsSumExp( x, _expX );
                    time_grad_wl += seconds() - time_used;
                    //             UpdatePotentialGrid(x);
                }// end cg
            }// end traversal nets

#if 0
            int merge_step = 1;
            if(ite % merge_step == 0){//
                double SL_step = 0;
                AdjustForce(n,x,SL_grad_f);
                LineSearch(n, x, SL_grad_f, SL_step);
                for(unsigned int i = 0; i < m_pDB->m_modules.size();i++){
                    x[2*i]   += SL_step*SL_grad_f[2*i];
                    x[2*i + 1] += SL_step*SL_grad_f[2*i + 1];
                    lastPosition[2*i] = x[2*i];
                    lastPosition[2*i+1] = x[2*i+1];
                }
                BoundX( n, x, x_l, x_u );
                UpdateExpValueForEachCell( n, x, _expX, _alpha );
                UpdateExpValueForEachPin( n, x, _expPins, _alpha );
                UpdateNetsSumExp( x, _expX );
                UpdateBlockPosition(x);
                UpdateDensityGrid(n,x);
                UpdatePotentialGrid(x);

            }
#endif
            printf("*");
        } // end while

        totalIte += innerIte;

        UpdateDensityGrid(n,x);
        UpdateBlockPosition(x);

        double nnb_real = GetNonZeroDensityGridPercent();
        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();
        if( param.bShow )
        {
            sprintf( filename, "SLfig%d-%d.plt", currentLevel, ite+1 );
            m_pDB->OutputGnuplotFigure( filename, false, false );
        }

        if(param.bShow){
            printf( "\n%d-%2d HPWL= %.0f  maxDen= %.2f OD= %.4f ODLB=%.4f OP=%.4f  LTime= %.1fm  objValue=%.0f, WireW= %.0f \n",
                    currentLevel, ite+1, m_pDB->CalcHPWL(),
                    maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                    double(seconds()-time_start)/60.0,
                    obj_value,
                    _weightWire
                    );
        }

        if (ite >= 2 && m_lookAheadLegalization && over < target_density + 0.1 )//2
        {
            //printf("in the loop,over < target_density + 0.1\n");
            UpdateBlockPosition(x);
            double hpwl = m_pDB->CalcHPWL();
            if (hpwl > bestLegalWL)
            {
                printf("Stop. Good enough.\n");
                break;
            }

            lookAheadLegalCount++;
            double oldWL = hpwl;
            CTetrisLegal legal(*m_pDB);

            double scale = 0.85;
            if (givenTargetUtil < 1.0 && givenTargetUtil > 0)
            {
                scale = 0.9;
            }

            double legalStart = seconds();
            bool   bLegal = legal.Solve( givenTargetUtil, false, false, scale);
            double legalTime = seconds() - legalStart;

            totalLegalTime += legalTime;

            if (param.bShow)
            {
                printf("LAL Time: %.2f\n", legalTime);
            }

            if (bLegal)
            {
                double WL = m_pDB->GetHPWLdensity(givenTargetUtil);
                if (param.bShow)
                {
                    m_pDB->ShowDensityInfo();
                }

                if (WL < bestLegalWL)
                {
                    LALnoGoodCount = 0;
                    if (param.bShow)
                    {
                        printf("SAVE BEST! (HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n",
                               m_pDB->GetHPWLp2p(), WL, (WL - oldWL) / oldWL * 100);
                    }
                    bestLegalWL = WL;
                    hasBestLegalSol = true;

                    for (int i = 0; i < (int) m_pDB->m_modules.size(); i++)
                    {
                        xBest[2 * i]     = m_pDB->m_modules[i].m_cx;
                        xBest[2 * i + 1] = m_pDB->m_modules[i].m_cy;
                        //                        x[2 * i]     = m_pDB->m_modules[i].m_cx;
                        //                        x[2 * i + 1] = m_pDB->m_modules[i].m_cy;

                    }
                }
                else {
                    if( param.bShow )
                        printf( "(HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n", m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    if( (WL-oldWL)/oldWL < 0.075 ) // 0.01
                    {
                        if( param.bShow )
                            printf( "Stop. Good enough\n" );
                        break;
                    }
                    LALnoGoodCount++;
                    if( LALnoGoodCount >= 2 )
                        break;
                }
            }
        } //end of m_lookAheadLegalization
        if (ite >= 3)//2  hdgao changed
        {
            if ( startDecreasing && over < target_density)// + 0.01
            {
                printf("Meet constraint!\n");
                break;
            }


            if (ite > 3 && totalOverPotential > lastTotalOverPotential && totalOverPotential < 1.4)
            {
                printf("Cannot further reduce!\n");
                break;
            }
        }


        _weightWire /= m_incFactor;   //  _weightDensity *=2;
        //       _weightDensity *=2;
        lastTotalOverPotential = totalOverPotential;

    }// end of "for (ite = 0; ite < maxIte; ite++)"

    if (hasBestLegalSol)
    {
        memcpy(x, xBest, sizeof(double) * n);
    }

    UpdateBlockPosition(x);

    if( lookAheadLegalCount > 0 && param.bShow )
    {
        printf( "LAL: Total Count: %d\n", lookAheadLegalCount );
        printf( "LAL: Total CPU: %.2f\n", totalLegalTime );
        sprintf( filename, "util-global.dat" );
        CPlaceBin placeBin( *m_pDB );
        placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
        placeBin.OutputBinUtil( filename );
    }

    delete [] lastPosition;


    delete [] SL_grad_f;
    return hasBestLegalSol;
}

bool   MyNLP::PLSolve_smooth(double wWire, double target_density, int currentLevel,threadpool_t* pool)
{
    double givenTargetUtil = m_targetUtil;    // m_targetUtil = -1
    m_currentStep = param.step;

    m_targetUtil += 0.05;//
    if (m_targetUtil > 1.0)
        m_targetUtil = 1.0;

    double time_start = seconds();
    char filename[100];    // for gnuplot

    size_t size_module_2 = 2 * m_pDB->m_modules.size();  //for the total num of modules
    // calculate the goal of the utilization
    double designUtil = m_pDB->m_totalModuleArea / m_pDB->m_totalFreeSpace;

    if (param.bShow)
        printf( "hdgao INFO: Design utilization: %f\n", designUtil);

    if (m_targetUtil > 0)
    {  // has give a utilization
        double lowest = designUtil + 0.05;
        if (m_targetUtil < lowest)
        {
            if (param.bShow)
            {
                printf("hdgao WARNING: Target utilization (%f) is too low\n", m_targetUtil);
                printf("hdgao          Set target utilization to %f \n", lowest);
            }
            m_targetUtil = lowest;
        }
    } else { // no given utilization
        printf("hdgao WARNING: No given target utilization. Distribute blocks evenly.\n");
        m_targetUtil = designUtil + 0.05;
        if (m_targetUtil > 1.0)
            m_targetUtil = 1.0;
    }

    if (param.bShow){
        printf("hdgao DBIN: param.run_cycle   : %f\n", param.run_cycle);
        printf("hdgao DBIN: Target utilization: %f\n", m_targetUtil);
    }




    double* lastPosition = new double[size_module_2];
    memset(lastPosition, 0, sizeof(double) * size_module_2);


    //calculate the subgrad of f and the direction and allocate memeory
    double* p_grad_f      = new double[ size_module_2 ];
    double* last_p_grad_f = new double[ size_module_2 ];
    memset(p_grad_f,     0,sizeof(double)*size_module_2);
    memset(last_p_grad_f,0,sizeof(double)*size_module_2);

    double* SL_grad_f = new double [size_module_2];
    memset( SL_grad_f, 0, sizeof(double)*size_module_2 );

    CreatePotentialGrid();   // create potential grid according to "m_potentialGridSize"
    int densityGridSize = 10;	// 1% chip area
    //int densityGridSize = m_potentialGridSize / 3;	    // not good in big3
    CreateDensityGrid( densityGridSize );	// real density: use 1% area
    //flpeng add so far the x has been calculated by qsolve or the last uncoarsening iteration
    UpdateDensityGridSpace( size_module_2, x );
    UpdatePotentialGridBase( x );		// init exp potential for each bin, also update ExpBin


#if 1
    // gaussian smoothing for base potential
    GaussianSmooth smooth;
    int r = m_smoothR;
    smooth.Gaussian2D( r, 6*r+1 );
    smooth.Smooth( m_basePotential );
    m_basePotentialOri = m_basePotential;
    sprintf( filename, "gbase%d.dat", currentLevel );
    OutputPotentialGrid( filename );
#endif

    // TEST
    if( m_smoothDelta == 1 )
    {
        if( param.bShow )
        {
            sprintf( filename, "gbase%d-more.dat", currentLevel );
            printf( "generate %s...\n", filename );
            fflush( stdout );
        }

        vector< vector< double > > moreSmooth = m_basePotential;
        r = m_smoothR * 6;
        int kernel_size = 5*r;
        if( kernel_size % 2 == 0 )
            kernel_size++;
        smooth.Gaussian2D( r, kernel_size );
        smooth.Smooth( moreSmooth );

        if( param.bShow )
        {
            swap( moreSmooth, m_basePotential );
            OutputPotentialGrid( filename );
            swap( moreSmooth, m_basePotential );
        }

        // merge base and moreSmooth
        double binArea = m_potentialGridWidth * m_potentialGridHeight;
        double halfBinArea = binArea / 2;
        int changeCount = 0;
        for( unsigned int i=0; i<moreSmooth.size(); i++ )
        {
            for( unsigned int j=0; j<moreSmooth[i].size(); j++ )
            {
                double free = binArea - m_basePotential[i][j];
                if( free < 1e-4 )	// no space
                {
                    if( moreSmooth[i][j] > halfBinArea )
                    {
                        m_basePotential[i][j] += moreSmooth[i][j] - halfBinArea;
                        changeCount++;
                    }
                }
            }
        }

        if( param.bShow )
        {
            printf( "change %d\n", changeCount );
            sprintf( filename, "gbase%d-more-merge.dat", currentLevel );
            OutputPotentialGrid( filename );
        }
    }


    if( m_smoothDelta > 1.0 )
        SmoothPotentialBase( double(m_smoothDelta) );   // also update ExpBin

    UpdateExpBinPotential( m_targetUtil );

    if( param.bShow )
    {
        sprintf( filename, "base%d.dat", currentLevel );
        OutputPotentialGrid( filename );
    }

    //   double objValue;

    assert( m_targetUtil > 0 );

    // wirelength
    UpdateExpValueForEachCell( size_module_2, x, _expX, _alpha );
    UpdateExpValueForEachPin( size_module_2, x, _expPins, _alpha );
    UpdateNetsSumExp( x, _expX );
    totalWL = GetWL( size_module_2, x, _expX, _alpha );

    // density
    UpdatePotentialGrid( x );
    UpdateDensityGrid( size_module_2, x );
    density = GetDensityPanelty();


    // 2006-02-22 weight (APlace ICCAD05)
    _weightWire = 1.0;
    eval_grad_f( size_module_2, x, _expX, true, p_grad_f );
    double totalWireGradient = 0;
    double totalPotentialGradient = 0;

    // TODO: truncation?
    AdjustForce( size_module_2, x, grad_wire, grad_potential );

    for( size_t i=0; i<size_module_2; i++ )
    {
        totalWireGradient      += fabs( grad_wire[i] );
        totalPotentialGradient += fabs( grad_potential[i] );
    }


    _weightDensity = 1.0;//1.02
    _weightWire = wWire *  totalPotentialGradient/totalWireGradient;//*param.sg_namdaIncFactor;//* param.sg_namdaIncFactor ;

    if(currentLevel != 1)
        _weightWire = _weightWire*1;

    int      maxIte = 50;   // max iterator
    //    double   beta;   // for the CG method
    double obj_value;
    eval_f(size_module_2,x,_expX,true,obj_value);
    eval_grad_f(size_module_2,x,_expX,true,p_grad_f);

    double nnbReal = GetNonZeroDensityGridPercent();
    UpdateDensityGrid( size_module_2, x );
    double maxDen = GetMaxDensity();
    double totalOverDen = GetTotalOverDensity();
    double totalOverDenLB = GetTotalOverDensityLB();
    double totalOverPotential = GetTotalOverPotential();

    if (param.bShow)
    {
        printf(" %d-%2d HPWL= %.0f\tDen= %.2f %.2f %.2f %.2f NNB= %.2f Dcost= %4.1f%%  WireW= %.0f",
               currentLevel, m_ite, m_pDB->CalcHPWL(), maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
               nnbReal, density * _weightDensity / obj_value * 100.0, _weightWire
               );
    } else {
        printf( " %d-%2d HPWL= %.0f \t", currentLevel, m_ite, m_pDB->CalcHPWL() );
    }
    fflush (stdout);

    if( param.bShow )
    {
        sprintf( filename, "PLfig%d-%d.plt", currentLevel, m_ite );
        m_pDB->OutputGnuplotFigure( filename, false, false );
    }


    double lastTotalOver = 0.0;
    double lastTotalOverPotential = DBL_MAX;
    double over = totalOverDen;
    //    int    totalIte = 0;
    bool   hasBestLegalSol = false;
    double bestLegalWL = DBL_MAX;
    int    lookAheadLegalCount = 0;
    double totalLegalTime = 0.0;

    bool   startDecreasing = (currentLevel == 1) ? false : true;

    int checkStep = 5;

    int LALnoGoodCount = 0;
    int totalIte = 0;

    int bucket_size = 10;
    double pick_rate = 1;//1: train all nets
    vector< vector<int> > temp_nets_cluster;
    temp_nets_cluster.resize(bucket_size * bucket_size);//temp_nets_cluster[0][0:99]
    vector< vector< vector<int> > > bucket;
    //test();
    if(1 && currentLevel != 1){
#if 0
        NCA_deeper();
        nets_cluster_modules_bounds.resize( nets_cluster.size() );
        for(int i = 0;i < nets_cluster.size();i++){
            getBlackModule(i);
        }
#endif

        double bucket_width  = (m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left)/bucket_size;
        double bucket_height = (m_pDB->m_coreRgn.top - m_pDB->m_coreRgn.bottom)/bucket_size;



        bucket.resize(bucket_size);
        for(int i = 0;i < bucket_size;i++)
            bucket[i].resize(bucket_size);

        //        for(int i = 0;i < bucket_size;i++)
        //            for(int j = 0;j < bucket_size;j++)
        //                bucket[i][j].clear();

        for(unsigned int i = 0;i < m_pDB->m_nets.size();i++){
            double c_x,c_y;
            get_center_position(x,c_x,c_y,i,0);//0: ortho_center,1: mean_center
            unsigned int n_x = (c_x - m_pDB->m_coreRgn.left)/(bucket_width );
            unsigned int n_y = (c_y - m_pDB->m_coreRgn.bottom)/(bucket_height);
            bucket[n_x][n_y].push_back(i);
        }
        // record all nets

        for(int i = 0;i < bucket_size;i++)
            for(int j = 0;j < bucket_size;j++){
                int num_pick_net = pick_rate*bucket[i][j].size();
                if(bucket[i][j].size() == 0) continue;
                if(num_pick_net == 0) num_pick_net = 1;
                vector<int> rand_vec = randVector(bucket[i][j].size());
                for(int k = 0;k < num_pick_net;k++)
                    temp_nets_cluster[i*bucket_size + j].push_back(bucket[i][j][ rand_vec[k] ]);
            }// as optimizer goes, nets may move out their original bucket


        if(param.bShow){
            printf("\n");
            int t = 0;
            for(int i = 0;i < bucket_size;i++){
                for(int j = 0;j < bucket_size;j++){
                    printf("%4d\t",bucket[i][j].size());// print bucket size
                    t += bucket[i][j].size();
                }
                printf("\n");
            }
            printf("net size  %d\n",t);
        }
        //        getchar();
    }
    else{//currentLevel = 1
        //whole net

        nets_cluster.resize(1);
        for(size_t i = 0;i < m_pDB->m_nets.size();i++)
            nets_cluster[0].push_back(i);
        nets_cluster_modules_bounds.resize( nets_cluster.size() );
        getBlackModule(0);
    }
    //   int cycle_limit = (currentLevel == 1) ? 100000:param.run_cycle;
    //    save_x_threads.resize(THREAD_size);// to record every results of multi-threads cg.
    double* x_bak = new double[size_module_2];
    bool is_2D = false;

    /*
    m_currentStep = param.step;//param.sg_stepsize
    if(currentLevel != 1){
        m_currentStep = 0.05;//param.run_cycle;
        PL_step = param.run_cycle;
        parallel_flag   = true;
       // param.run_cycle = PL_step;
        printf("PL_step = %f\n",PL_step);
      //  getchar();
    }
*/

    //    MyNLP *thread_mynlp[THREAD_size];
    //    args_class *ac[THREAD_size];
    vector<args_class*> ac;
    ac.reserve(THREAD_size);
    size_t avg_task = bucket_size/THREAD_size;
    if(currentLevel != 1){
        for(int i = 0; i < THREAD_size;i++){

            MyNLP* temp = new MyNLP;
            temp->MyNLP_copy(this);

            vector<int> specified_nets_cluster_Id;
            for(size_t k = 0;k < avg_task;k++)// just for case: Thread_size = 1
                specified_nets_cluster_Id.push_back(i*avg_task + k);
            ac[i] = new args_class(temp,//thread_mynlp[i]
                                   specified_nets_cluster_Id,size_module_2,
                                   nets_cluster,
                                   nets_black_modules,
                                   nets_cluster_modules_bounds,
                                   is_2D);
        }
        printf("ac done\n");
    }
    /*
        {
            double pick_rate = 0.01;
            vector< vector<int> > temp_nets_cluster;
            nets_cluster.clear();
            nets_black_modules.clear();
            nets_cluster_modules_bounds.clear();

            temp_nets_cluster.resize(bucket_size * bucket_size);//temp_nets_cluster[0][0:99]

            for(int i = 0;i < bucket_size;i++)
                for(int j = 0;j < bucket_size;j++){
                    int num_pick_net = pick_rate*bucket[i][j].size();
                    if(bucket[i][j].size() == 0) continue;
                    if(num_pick_net == 0) num_pick_net = 1;
                    vector<int> rand_vec = randVector(bucket[i][j].size());
                    for(int k = 0;k < num_pick_net;k++)
                        temp_nets_cluster[i*bucket_size + j].push_back(bucket[i][j][ rand_vec[k] ]);
                }
            vector<int> rand_vec = randVector(bucket_size*bucket_size);

            for(int m = 0; m < bucket_size;m++){// 100 buckets in total.
                vector<int> temp;
                for(int i = 0;i < bucket_size;i++)// each nets_cluster contains 10 buckets
                    temp.insert(temp.end(),temp_nets_cluster[ rand_vec[ m*bucket_size + i] ].begin(),
                                           temp_nets_cluster[ rand_vec[ m*bucket_size + i] ].end() );
                nets_cluster.push_back(temp);// nets_cluster[0][...]
            }
            nets_black_modules.resize(nets_cluster.size());
            nets_cluster_modules_bounds.resize(nets_cluster.size());
            for(unsigned int i = 0; i < nets_cluster.size();i++)
                getBlackModule(i);

            vector<int> specified_nets_cluster_Id;
            specified_nets_cluster_Id.push_back(0);
            args_class* ac = new args_class(thread_mynlp[0],specified_nets_cluster_Id,n,
                                            nets_cluster,nets_black_modules,nets_cluster_modules_bounds,
                                            false);
        }
 */

//    pthread_mutex_init(&lock, nullptr);
    job_done = 0;
#if 0
    test();
#endif

    double* x_bak_bak = new double[size_module_2];
    for(int i = 0; i < size_module_2;i++){
        x_bak_bak[i] = x[i];
    }

    printf("\n======================================================================================\n");
    printf("cg begin\n");

    for(int ite = 0;ite < maxIte;ite++){
        if( _weightWire == 0) break;
        int innerIte = 0;
        double old_obj = DBL_MAX;
        double last_obj_value = DBL_MAX;

        for(int i = 0; i < size_module_2;i++){// to update potential more quickly
            x_bak[i] = x[i];
        }

        while( (currentLevel == 1) || (currentLevel != 1 && innerIte <= 1) ){
            innerIte++;


            if( innerIte % checkStep == 0){
                old_obj = last_obj_value;
                eval_f(size_module_2,x,_expX,true,obj_value);
                last_obj_value = obj_value;
            }
            if( innerIte % checkStep == 0){
                printf(".");
                fflush( stdout );
                if( innerIte % (2*checkStep) == 0){
                    UpdateBlockPosition( x );
                    if( m_pDB->CalcHPWL() > bestLegalWL ){
                        printf("best legal WL\n");
                        break;
                    }
                }

                UpdateDensityGrid(size_module_2,x);
                totalOverDen = GetTotalOverDensity();
                totalOverDenLB = GetTotalOverDensityLB();
                totalOverPotential = GetTotalOverPotential();

                lastTotalOver = over;
                over = min( totalOverPotential, totalOverDen ); // TEST

                if( !startDecreasing
                        && over < lastTotalOver
                        && ite >= 1
                        && innerIte >= 6)
                {
                    printf( ">>" );
                    fflush( stdout );
                    startDecreasing = true;
                }

                // 2005-03-11: meet the constraint
                if( startDecreasing && over < target_density )
                    break;
                if( obj_value >= param.precision * old_obj){
                    break;
                }
            }

            /*  pick random net */
            if(currentLevel != 1 && innerIte == 1){


                nets_cluster.clear();
                nets_black_modules.clear();
                nets_cluster_modules_bounds.clear();


                vector<int> rand_vec = randVector(bucket_size*bucket_size);

                for(int m = 0; m < bucket_size;m++){// 100 buckets in total.
                    vector<int> temp;
                    for(int i = 0;i < bucket_size;i++)// each nets_cluster contains 10 buckets
                        temp.insert(temp.end(),temp_nets_cluster[ rand_vec[ m*bucket_size + i] ].begin(),
                                temp_nets_cluster[ rand_vec[ m*bucket_size + i] ].end() );
                    nets_cluster.push_back(temp);// nets_cluster[0][...]
                }
                // nets_black_modules.resize(nets_cluster.size());
                nets_cluster_modules_bounds.resize(nets_cluster.size());
                for(unsigned int i = 0; i < nets_cluster.size();i++)
                    getBlackModule(i);
            }

            //  update modules positions
            if(currentLevel == 1){
                //              CG_mynlp(n,x_bak,is_2D);
                CG(this,size_module_2,x_bak,is_2D);
            }
            else{// parallel
                job_done = 0;
                done_flag = false;

                //update cluster
                for(unsigned int i_thread = 0;i_thread < THREAD_size;i_thread++){
                    ac[i_thread]->args_mynlp->nets_cluster.clear();
                    ac[i_thread]->args_mynlp->nets_black_modules.clear();
                    ac[i_thread]->args_mynlp->nets_cluster_modules_bounds.clear();

                    memcpy(ac[i_thread]->args_mynlp->x,x,sizeof(double)*size_module_2);
                    //   for(int i = 0;i < n;i++) ac[i_thread]->args_mynlp->x[i] = x[i];
                    ac[i_thread]->args_mynlp->nets_cluster.resize(nets_cluster.size());
                    ac[i_thread]->args_mynlp->nets_black_modules.resize(nets_black_modules.size());
                    ac[i_thread]->args_mynlp->nets_cluster_modules_bounds.resize(nets_cluster_modules_bounds.size());

                    ac[i_thread]->args_mynlp->_weightDensity = _weightDensity;
                    ac[i_thread]->args_mynlp->_weightWire = _weightWire;

                    for(size_t i = 0;i < nets_cluster.size();i++)
                        ac[i_thread]->args_mynlp->nets_cluster[i].assign(nets_cluster[i].begin(),
                                                                         nets_cluster[i].end());


                    for(size_t i = 0;i < nets_black_modules.size();i++)
                        ac[i_thread]->args_mynlp->nets_black_modules[i].assign(nets_black_modules[i].begin(),
                                                                               nets_black_modules[i].end());

                    for(size_t i = 0;i < nets_cluster_modules_bounds.size();i++)
                        ac[i_thread]->args_mynlp->nets_cluster_modules_bounds[i].assign(nets_cluster_modules_bounds[i].begin(),
                                                                                        nets_cluster_modules_bounds[i].end() );

                    ac[i_thread]->thread_number = i_thread;
                    assert( threadpool_add(pool,MyNLP::task_GD,(void*)(ac[i_thread]),0) == 0 );
                }

                while(1){
                    if(job_done == THREAD_size || save_x_threads.size() == THREAD_size){
                        break;
                    }
                    pthread_mutex_unlock(&lock);
                }
            }
            // merge
#if 1
            if(currentLevel != 1){
                int merge_step = 1;

                if(ite % merge_step == 0){//
                    for(size_t j = 0;j < size_module_2/2;j++){
                        double x_sum = 0;
                        double y_sum = 0;
                        for(size_t i = 0;i < THREAD_size;i++){
                            x_sum += save_x_threads[i][ 2*j ];
                            y_sum += save_x_threads[i][ 2*j+1 ];
                        }
                        x[2*j]   = x_sum/THREAD_size;
                        x[2*j+1] = y_sum/THREAD_size;
                    }
                    //print
#if 1
                    ofstream x_file;
                    x_file.open("x_record.txt",ostream::app);
                    x_file << "loop: " << ite << endl;
                    for(int i =0;i < size_module_2;i++){
                        for(int j = 0;j < THREAD_size;j++){
                            x_file << save_x_threads[j][i] <<"\t";
                        }
                        x_file << x[i] << "\t";
                        x_file << x_bak_bak[i] << "\t";
                        x_file << fabs(x[i] - x_bak_bak[i]) << endl;
                    }
                    x_file.close();
#endif

                    BoundX( size_module_2, x, x_l, x_u );
                    UpdateExpValueForEachCell( size_module_2, x, _expX, _alpha );
                    UpdateExpValueForEachPin( size_module_2, x, _expPins, _alpha );
                    UpdateNetsSumExp( x, _expX );
                    UpdatePotentialGrid(x);

                    //                    UpdateBlockPosition(x);
                    //                  UpdateDensityGrid(n,x);

                    save_x_threads.clear();
                    /*
                        ac[i_thread]->args_mynlp->BoundX(n, ac[i_thread]->args_mynlp->x,
                                                         ac[i_thread]->args_mynlp->x_l, ac[i_thread]->args_mynlp->x_u);
                        ac[i_thread]->args_mynlp->UpdateExpValueForEachCell(n, ac[i_thread]->args_mynlp->x,
                                                                            ac[i_thread]->args_mynlp->_expX, ac[i_thread]->args_mynlp->_alpha);
                        ac[i_thread]->args_mynlp->UpdateExpValueForEachPin(n, ac[i_thread]->args_mynlp->x,
                                                                           ac[i_thread]->args_mynlp->_expX, ac[i_thread]->args_mynlp->_alpha);
                        ac[i_thread]->args_mynlp->UpdateNetsSumExp(ac[i_thread]->args_mynlp->x,
                                                                   ac[i_thread]->args_mynlp->_expX);
                        ac[i_thread]->args_mynlp->UpdateBlockPosition(ac[i_thread]->args_mynlp->x);
                        ac[i_thread]->args_mynlp->UpdateDensityGrid(n,ac[i_thread]->args_mynlp->x);
                        ac[i_thread]->args_mynlp->UpdatePotentialGrid(ac[i_thread]->args_mynlp->x);
                        */

                }
            }
#endif
            if(currentLevel != 1 ) printf("*");
            fflush(stdout);
        } // end while

        totalIte += innerIte;

        UpdateDensityGrid(size_module_2,x);
        UpdateBlockPosition(x);

        //        double nnb_real = GetNonZeroDensityGridPercent();
        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();
        over = min( totalOverPotential, totalOverDen ); // TEST
        if( param.bShow )
        {
            if(ite+1 < 10){
                sprintf( filename, "PLfig%d-0%d.plt", currentLevel, ite+1 );
                m_pDB->OutputGnuplotFigure( filename, false, false );
            }
            else{
                sprintf( filename, "PLfig%d-%d.plt", currentLevel, ite+1 );
                m_pDB->OutputGnuplotFigure( filename, false, false );
            }
        }

        if(param.bShow){
            printf( "\n%d-%2d HPWL= %.0f  maxDen= %.2f OD= %.4f ODLB=%.4f OP=%.4f  LTime= %.1fm  objValue=%.0f, WireW= %.0f \n",
                    currentLevel, ite+1, m_pDB->CalcHPWL(),
                    maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                    double(seconds()-time_start)/60.0,
                    obj_value,
                    _weightWire
                    );
        }

        if (1 && ite >= 2 && m_lookAheadLegalization && over < target_density + 0.10 )//2
        {
            //printf("in the loop,over < target_density + 0.1\n");
            UpdateBlockPosition(x);
            double hpwl = m_pDB->CalcHPWL();
            if (hpwl > bestLegalWL)
            {
                printf("Stop. Good enough 1.\n");
                break;
            }

            lookAheadLegalCount++;
            double oldWL = hpwl;
            CTetrisLegal legal(*m_pDB);

            double scale = 0.85;
            if (givenTargetUtil < 1.0 && givenTargetUtil > 0)
            {
                scale = 0.9;
            }

            double legalStart = seconds();
            bool   bLegal = legal.Solve( givenTargetUtil, false, false, scale);
            double legalTime = seconds() - legalStart;

            totalLegalTime += legalTime;

            if (param.bShow)
            {
                printf("LAL Time: %.2f\n", legalTime);
            }

            if (bLegal)
            {
                double WL = m_pDB->GetHPWLdensity(givenTargetUtil);
                if (param.bShow)
                {
                    m_pDB->ShowDensityInfo();
                }

                if (WL < bestLegalWL)
                {
                    LALnoGoodCount = 0;
                    if (param.bShow)
                    {
                        printf("SAVE BEST! (HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n",
                               m_pDB->GetHPWLp2p(), WL, (WL - oldWL) / oldWL * 100);
                    }
                    bestLegalWL = WL;
                    hasBestLegalSol = true;

                    for (size_t i = 0; i <  m_pDB->m_modules.size(); i++)
                    {
                        xBest[2 * i]     = m_pDB->m_modules[i].m_cx;
                        xBest[2 * i + 1] = m_pDB->m_modules[i].m_cy;
                        //                        x[2 * i]     = m_pDB->m_modules[i].m_cx;
                        //                        x[2 * i + 1] = m_pDB->m_modules[i].m_cy;

                    }
                }
                else {
                    if( param.bShow )
                        printf( "(HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n", m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    if( (WL-oldWL)/oldWL < 0.075 ) // 0.01
                    {
                        if( param.bShow )
                            printf( "Stop. Good enough 2\n" );
                        break;
                    }
                    LALnoGoodCount++;
                    if( LALnoGoodCount >= 2 ){
                        printf("LALnoGoodCount >= 2\n");
                        break;
                    }
                }
            }
        } //end of m_lookAheadLegalization
        if (ite >= 2)//2  hdgao changed
        {
            if ( startDecreasing && over < target_density)// + 0.01
            {
                printf("Meet constraint!\n");
                break;
            }


            if (ite > 3 && totalOverPotential > lastTotalOverPotential &&
                    totalOverPotential < 1.4)
            {
                printf("Cannot further reduce!\n");
                break;
            }
        }


        _weightWire /= m_incFactor;   //  _weightDensity *=2;
        //       _weightDensity *=2;
        lastTotalOverPotential = totalOverPotential;

    }// end of "for (ite = 0; ite < maxIte; ite++)"

    if (hasBestLegalSol)
    {
        memcpy(x, xBest, sizeof(double) * size_module_2);
    }

    UpdateBlockPosition(x);

    if( lookAheadLegalCount > 0 && param.bShow )
    {
        printf( "LAL: Total Count: %d\n", lookAheadLegalCount );
        printf( "LAL: Total CPU: %.2f\n", totalLegalTime );
        sprintf( filename, "util-global.dat" );
        CPlaceBin placeBin( *m_pDB );
        placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
        placeBin.OutputBinUtil( filename );
    }

    delete [] lastPosition;

    delete [] SL_grad_f;
    return hasBestLegalSol;
}
bool   MyNLP::ParallelSolve(double wWire, double target_density, int currentLevel,threadpool_t** pool)
{
    double givenTargetUtil = m_targetUtil;    // m_targetUtil = -1
    m_currentStep = param.step;

    m_targetUtil += 0.05;//
    if (m_targetUtil > 1.0)
        m_targetUtil = 1.0;

    double time_start = seconds();
    char filename[100];    // for gnuplot

    size_t size_module_2 = 2 * m_pDB->m_modules.size();  //for the total num of modules
    // calculate the goal of the utilization
    double designUtil = m_pDB->m_totalModuleArea / m_pDB->m_totalFreeSpace;

    if (param.bShow)
        printf( "hdgao INFO: Design utilization: %f\n", designUtil);

    if (m_targetUtil > 0)
    {  // has give a utilization
        double lowest = designUtil + 0.05;
        if (m_targetUtil < lowest)
        {
            if (param.bShow)
            {
                printf("hdgao WARNING: Target utilization (%f) is too low\n", m_targetUtil);
                printf("hdgao          Set target utilization to %f \n", lowest);
            }
            m_targetUtil = lowest;
        }
    } else { // no given utilization
        printf("hdgao WARNING: No given target utilization. Distribute blocks evenly.\n");
        m_targetUtil = designUtil + 0.05;
        if (m_targetUtil > 1.0)
            m_targetUtil = 1.0;
    }

    if (param.bShow){
        printf("hdgao DBIN: param.run_cycle   : %f\n", param.run_cycle);
        printf("hdgao DBIN: Target utilization: %f\n", m_targetUtil);
    }




    double* lastPosition = new double[size_module_2];
    memset(lastPosition, 0, sizeof(double) * size_module_2);


    //calculate the subgrad of f and the direction and allocate memeory
    double* p_grad_f      = new double[ size_module_2 ];
    double* last_p_grad_f = new double[ size_module_2 ];
    memset(p_grad_f,     0,sizeof(double)*size_module_2);
    memset(last_p_grad_f,0,sizeof(double)*size_module_2);


    CreatePotentialGrid();   // create potential grid according to "m_potentialGridSize"
    int densityGridSize = 10;	// 1% chip area
    //int densityGridSize = m_potentialGridSize / 3;	    // not good in big3
    CreateDensityGrid( densityGridSize );	// real density: use 1% area
    //flpeng add so far the x has been calculated by qsolve or the last uncoarsening iteration
    UpdateDensityGridSpace( size_module_2, x );
    UpdatePotentialGridBase( x );		// init exp potential for each bin, also update ExpBin


#if 1
    // gaussian smoothing for base potential
    GaussianSmooth smooth;
    int r = m_smoothR;
    smooth.Gaussian2D( r, 6*r+1 );
    smooth.Smooth( m_basePotential );
    m_basePotentialOri = m_basePotential;
    sprintf( filename, "gbase%d.dat", currentLevel );
    OutputPotentialGrid( filename );
#endif

    // TEST
    if( m_smoothDelta == 1 )
    {
        if( param.bShow )
        {
            sprintf( filename, "gbase%d-more.dat", currentLevel );
            printf( "generate %s...\n", filename );
            fflush( stdout );
        }

        vector< vector< double > > moreSmooth = m_basePotential;
        r = m_smoothR * 6;
        int kernel_size = 5*r;
        if( kernel_size % 2 == 0 )
            kernel_size++;
        smooth.Gaussian2D( r, kernel_size );
        smooth.Smooth( moreSmooth );

        if( param.bShow )
        {
            swap( moreSmooth, m_basePotential );
            OutputPotentialGrid( filename );
            swap( moreSmooth, m_basePotential );
        }

        // merge base and moreSmooth
        double binArea = m_potentialGridWidth * m_potentialGridHeight;
        double halfBinArea = binArea / 2;
        int changeCount = 0;
        for( unsigned int i=0; i<moreSmooth.size(); i++ )
        {
            for( unsigned int j=0; j<moreSmooth[i].size(); j++ )
            {
                double free = binArea - m_basePotential[i][j];
                if( free < 1e-4 )	// no space
                {
                    if( moreSmooth[i][j] > halfBinArea )
                    {
                        m_basePotential[i][j] += moreSmooth[i][j] - halfBinArea;
                        changeCount++;
                    }
                }
            }
        }

        if( param.bShow )
        {
            printf( "change %d\n", changeCount );
            sprintf( filename, "gbase%d-more-merge.dat", currentLevel );
            OutputPotentialGrid( filename );
        }
    }


    if( m_smoothDelta > 1.0 )
        SmoothPotentialBase( double(m_smoothDelta) );   // also update ExpBin

    UpdateExpBinPotential( m_targetUtil );

    if( param.bShow )
    {
        sprintf( filename, "base%d.dat", currentLevel );
        OutputPotentialGrid( filename );
    }

    //   double objValue;

    assert( m_targetUtil > 0 );

    // wirelength
    UpdateExpValueForEachCell( size_module_2, x, _expX, _alpha );
    UpdateExpValueForEachPin( size_module_2, x, _expPins, _alpha );
    UpdateNetsSumExp( x, _expX );
    totalWL = GetWL( size_module_2, x, _expX, _alpha );

    // density
    UpdatePotentialGrid( x );
    UpdateDensityGrid( size_module_2, x );
    density = GetDensityPanelty();

    // 2006-02-22 weight (APlace ICCAD05)
    _weightWire = 1.0;
    eval_grad_f( size_module_2, x, _expX, true, p_grad_f );
    double totalWireGradient = 0;
    double totalPotentialGradient = 0;

    // TODO: truncation?
    AdjustForce( size_module_2, x, grad_wire, grad_potential );

    for( size_t i=0; i<size_module_2; i++ )
    {
        totalWireGradient      += fabs( grad_wire[i] );
        totalPotentialGradient += fabs( grad_potential[i] );
    }


    _weightDensity = 1.0;//1.02
    _weightWire = wWire *  totalPotentialGradient/totalWireGradient;//*param.sg_namdaIncFactor;//* param.sg_namdaIncFactor ;

    double parallelScale = (param.isParallel)?param.loop*(param.Den_thread_size + param.WL_thread_size):1;
    //_weightDensity *= (parallelScale);
    //_weightWire *= parallelScale;

    if(currentLevel != 1)
        _weightWire = _weightWire*1;

    int      maxIte = 50;   // max iterator
    //    double   beta;   // for the CG method
    double obj_value;
    eval_f(size_module_2,x,_expX,true,obj_value);
    eval_grad_f(size_module_2,x,_expX,true,p_grad_f);

    double nnbReal = GetNonZeroDensityGridPercent();
    UpdateDensityGrid( size_module_2, x );
    double maxDen = GetMaxDensity();
    double totalOverDen = GetTotalOverDensity();
    double totalOverDenLB = GetTotalOverDensityLB();
    double totalOverPotential = GetTotalOverPotential();

    if (param.bShow)
    {
        printf(" %d-%2d HPWL= %.0f\tDen= %.2f %.2f %.2f %.2f NNB= %.2f Dcost= %4.1f%%  WireW= %.0f",
               currentLevel, m_ite, m_pDB->CalcHPWL(), maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
               nnbReal, density * _weightDensity / obj_value * 100.0, _weightWire
               );
    } else {
        printf( " %d-%2d HPWL= %.0f \t", currentLevel, m_ite, m_pDB->CalcHPWL() );
    }
    fflush (stdout);

    if( param.bShow )
    {
        sprintf( filename, "PLfig%d-%d.plt", currentLevel, m_ite );
        m_pDB->OutputGnuplotFigure( filename, false, false );
    }


    double lastTotalOver = 0.0;
    double lastTotalOverPotential = DBL_MAX;
    double over = totalOverDen;
    //    int    totalIte = 0;
    bool   hasBestLegalSol = false;
    double bestLegalWL = DBL_MAX;
    int    lookAheadLegalCount = 0;
    double totalLegalTime = 0.0;

    bool   startDecreasing = false;

    int checkStep = param.checkstep;

    int LALnoGoodCount = 0;
    int totalIte = 0;


    // data preparation
    bool isParallel = param.isParallel;
  //  param.Den_thread_size = 100;//100;//m_potentialGridSize*m_potentialGridSize;
   // param.WL_thread_size  = 100;//100;
    vector<args_class*> ac_Den;
    vector<args_class*> ac_WL;
    double bucket_width = 0,bucket_height = 0;
    size_t tot_tasks = param.Den_thread_size + param.WL_thread_size;
    size_t bin_size = sqrt(param.Den_thread_size);
    size_t net_size = m_pDB->m_nets.size()/param.WL_thread_size;
    if( isParallel ){
        ac_Den.resize(param.Den_thread_size);
        ac_WL.resize(param.WL_thread_size);
        bucket_width  = ( m_pDB->m_coreRgn.right - m_pDB->m_coreRgn.left)/bin_size;
        bucket_height = (m_pDB->m_coreRgn.top  - m_pDB->m_coreRgn.bottom)/bin_size;
        // for density
        nets_cluster_Den.resize(bin_size);
        for(size_t i = 0;i < bin_size;i++)
            nets_cluster_Den[i].resize(bin_size);
        for(unsigned int i = 0;i < m_pDB->m_nets.size();i++){
            double c_x,c_y;
            get_center_position(x,c_x,c_y,i,0);//0: ortho_center,1: mean_center
            unsigned int n_x = (c_x - m_pDB->m_coreRgn.left)/(bucket_width );
            unsigned int n_y = (c_y - m_pDB->m_coreRgn.bottom)/(bucket_height);
            nets_cluster_Den[n_x][n_y].push_back(i);
        }
        // for wirelength
        vector<int> rand_vec = randVector(m_pDB->m_nets.size());
        nets_cluster_WL.resize(param.WL_thread_size);
        vector<int>::iterator iter_temp = rand_vec.begin();
        for(size_t i = 0;i < param.WL_thread_size;i++){
            if(i != param.WL_thread_size - 1){
                nets_cluster_WL[i].assign( iter_temp + i*net_size,
                                           iter_temp + i*net_size + net_size);
            }
            else{// i = param.WL_thread_size - 1 = 3
                nets_cluster_WL[i].assign(iter_temp + i*net_size,rand_vec.end());
            }
        }
        for(size_t i = 0;i < bin_size;i++)
            for(size_t j = 0;j < bin_size; j++){// for density
                MyNLP* temp = new MyNLP( *(this->m_pDB) );
                {// for MyNLP
                    temp->m_potentialGridSize   = m_potentialGridSize;
                    temp->m_potentialGridHeight = m_potentialGridHeight;
                    temp->m_potentialGridWidth  = m_potentialGridWidth;
                    temp->_potentialRX = _potentialRX;
                    temp->_potentialRY = _potentialRY;
                    memcpy(temp->x_l,x_l,sizeof(double)*size_module_2);
                    memcpy(temp->x_u,x_u,sizeof(double)*size_module_2);
                    temp->xBak = new double[size_module_2];
                    memset(temp->xBak,0,sizeof(double)*size_module_2);
                    temp->_alpha = _alpha;
                    temp->_weightWire = _weightWire;
                    temp->_weightDensity = _weightDensity;

                    temp->m_gridPotential.resize(temp->m_potentialGridSize );
                    for(size_t i = 0;i < temp->m_potentialGridSize;i++)
                        temp->m_gridPotential[i].resize(temp->m_potentialGridSize );

                    temp->m_expBinPotential.resize(temp->m_potentialGridSize);
                    for(size_t i = 0;i < temp->m_potentialGridSize;i++)
                        temp->m_expBinPotential[i].assign(m_expBinPotential[i].begin(),
                                                          m_expBinPotential[i].end());

                }
                // for auxiliary parameters
                vector<int> specified_nets_cluster_Id;                
                specified_nets_cluster_Id.assign(nets_cluster_Den[i][j].begin(),
                                                 nets_cluster_Den[i][j].end());
                GetModulesByNets(specified_nets_cluster_Id,temp->modules);
                size_t i_thread = i*bin_size + j;
                ac_Den[i_thread] = new args_class(temp,
                                       specified_nets_cluster_Id,
                                       size_module_2);
                {                  
                    ac_Den[i_thread]->thread_number = i_thread;
                    //ac_Den[i_thread]->stepsize = 0.92;
                    ac_Den[i_thread]->maxloop  = param.loop;
                    ac_Den[i_thread]->args_mynlp->GetGrid(ac_Den[i_thread]->thread_number/bin_size,ac_Den[i_thread]->thread_number%bin_size,\
                                                          bucket_width,bucket_height,ac_Den[i_thread]->bin);
                }

        }
#if 0
        for(size_t i = 0;i < bin_size;i++)
            for(size_t j = 0;j < bin_size; j++){// for density
                MyNLP* temp = new MyNLP( *(this->m_pDB) );
                {// for MyNLP
                    temp->m_potentialGridSize   = m_potentialGridSize;
                    temp->m_potentialGridHeight = m_potentialGridHeight;
                    temp->m_potentialGridWidth  = m_potentialGridWidth;
                    temp->_potentialRX = _potentialRX;
                    temp->_potentialRY = _potentialRY;
                    memcpy(temp->x_l,x_l,sizeof(double)*size_module_2);
                    memcpy(temp->x_u,x_u,sizeof(double)*size_module_2);
                    //temp->xBak = new double[size_module_2];
                   // memset(temp->xBak,0,sizeof(double)*size_module_2);
                    temp->_alpha = _alpha;
                    temp->_weightWire = _weightWire;
                    temp->_weightDensity = _weightDensity;

                }
                // for auxiliary parameters
                vector<int> specified_nets_cluster_Id;
                specified_nets_cluster_Id.assign(nets_cluster_Den[i][j].begin(),
                                                 nets_cluster_Den[i][j].end());
                GetModulesByNets(specified_nets_cluster_Id,temp->modules);
                size_t i_thread = i*bin_size + j;
                ac_WL[i_thread] = new args_class(temp,
                                       specified_nets_cluster_Id,
                                       size_module_2);
                {
                    ac_WL[i_thread]->thread_number = i_thread + param.Den_thread_size;
                    ac_WL[i_thread]->maxloop  = param.loop;

                }

        }
#endif
#if 1
        for(size_t i = 0;i < param.WL_thread_size;i++){// for wirelength
            MyNLP* temp = new MyNLP( *(this->m_pDB) );
            {
                temp->m_potentialGridSize   = m_potentialGridSize;
                temp->m_potentialGridHeight = m_potentialGridHeight;
                temp->m_potentialGridWidth  = m_potentialGridWidth;
                temp->_potentialRX = _potentialRX;
                temp->_potentialRY = _potentialRY;
                memcpy(temp->x_l,x_l,sizeof(double)*size_module_2);
                memcpy(temp->x_u,x_u,sizeof(double)*size_module_2);
                //temp->xBak = new double[size_module_2];
                //memset(temp->xBak,0,size_module_2);
                temp->_alpha = _alpha;
                temp->_weightWire = _weightWire;
                temp->_weightDensity = _weightDensity;
                temp->m_gridPotential.resize(temp->m_potentialGridSize );
                for(size_t i = 0;i < temp->m_potentialGridSize;i++)
                    temp->m_gridPotential[i].resize(temp->m_potentialGridSize );
            }
            vector<int> specified_nets_cluster_Id;
            specified_nets_cluster_Id.assign(nets_cluster_WL[i].begin(),
                                             nets_cluster_WL[i].end());
            GetModulesByNets(specified_nets_cluster_Id,temp->modules);
            ac_WL[i] = new args_class(temp,
                                       specified_nets_cluster_Id,
                                       size_module_2);
            {                
                ac_WL[i]->thread_number = i + param.Den_thread_size;
                //ac_WL[i]->stepsize = 0.92;
                ac_WL[i]->maxloop = param.loop;
            }
        }
#endif

    }// end data preparation

    double* x_bak = new double[size_module_2];

    double* grad_f = new double [ 2 * m_pDB->m_modules.size() ];
    double* last_grad_f = new double [ 2 * m_pDB->m_modules.size() ]; // for computing CG-direction;
    memset( grad_f, 0, sizeof(double)*2*m_pDB->m_modules.size() );
    memset( last_grad_f, 0, sizeof(double)*2*m_pDB->m_modules.size() );




    pthread_mutex_init(&lock, nullptr);
#if 0
    if(currentLevel != 1){
        printf("\n");
        vector<double> grad_density,grad_wirelength;
        grad_density.resize(size_module_2,0);
        grad_wirelength.resize(size_module_2,0);
 /*        double sum1 = 0,sum2 = 0;
        double gradDensityX,gradDensityY,gradx,grady;
        int N = 0;
       for(size_t i = 0; i < size_module_2/2;i++){

            GetPotentialGrad( x, i, gradDensityX, gradDensityY );
            GetDensityGrad(i,ac_Den[0]->bin,gradx,grady);

            double delta1 = gradDensityX - gradx;
            double delta2 = gradDensityY - grady;
            sum1 += fabs(delta1);
            sum2 += fabs(delta2);
           if(delta1 > 0.001 || delta2 > 0.001){
                N++;
                printf("module i = %5d,gradDensityX = %.4f, gradDensityY = %.4f\n",
                       i,gradDensityX,gradDensityY);
                printf("module i = %5d,gradDensityX = %.4f, gradDensityY = %.4f\n",
                       i,gradx,grady);
                printf("delat1 = %.4f,delta2 = %.4f\n",delta1,delta2);

            }

        }
        printf("N = %d,sum1 = %.4f,sum2 = %.4f\n",N,sum1,sum2);
        getchar();
*/
        eval_grad_f( size_module_2, x, _expX, true, grad_f );
        DensityGrad(grad_density,ac_Den[0]->args_snc_Id,ac_Den[0]->bin);
        WireLengthGrad(grad_wirelength,ac_WL[0]->args_snc_Id);
        double sum_diff = 0;
        for(size_t i = 0;i < size_module_2;i++){
            double temp = grad_density[i]*1 +
                    grad_wirelength[i]*1;
            double delta =  grad_f[i] - temp;
            printf("gradx = %.4f, grad_potential = %.4f, grad_wire = %.4f\n",
                   grad_f[i],grad_potential[i]*_weightDensity,grad_wire[i]*_weightWire);
            printf("gradx = %.4f, grad_potential = %.4f, grad_wire = %.4f\ndelta = %.4f\n",
                   temp,grad_density[i],
                   grad_wirelength[i],
                   delta);
            sum_diff += delta;
       //     getchar();
        }
        printf("sum_diff = %4f\n",sum_diff);
        getchar();

}
#endif
    //for test
    vector<double> grad_density,grad_wirelength,grad;
    bool GD = param.gd;
    vector<int> modules;
    vector<int> whole_nets;
    double stepsize;
    vector<int> bin = {0,0,m_potentialGridSize-1,m_potentialGridSize-1};
if(1){
    grad_density.resize(size_module_2,0);
    grad_wirelength.resize(size_module_2,0);
    grad.resize(size_module_2,0);
    stepsize = 0.92;
    for(size_t i = 0;i < m_pDB->m_nets.size();i++)
        whole_nets.push_back(i);
    memcpy(xBak,x,sizeof(double)*size_module_2);
    if(0){
        eval_grad_f( size_module_2, x, _expX, true, grad_f );
        DensityGrad(grad_density,whole_nets,bin);
        WireLengthGrad(grad_wirelength,whole_nets);
        double sum_diff = 0;
        for(size_t i = 0;i < size_module_2;i++){
            double temp = grad_density[i]*1 +
                    grad_wirelength[i]*1;
            double delta =  grad_f[i] - temp;
            printf("gradx = %.4f, grad_potential = %.4f, grad_wire = %.4f\n",
                   grad_f[i],grad_potential[i]*_weightDensity,grad_wire[i]*_weightWire);
            printf("gradx = %.4f, grad_potential = %.4f, grad_wire = %.4f\ndelta = %.4f\n",
                   temp,grad_density[i],
                   grad_wirelength[i],
                   delta);
            sum_diff += delta;
       //     getchar();
        }
        printf("sum_diff = %4f\n",sum_diff);
        getchar();
    }

    GetModulesByNets(whole_nets,modules);
}
#if 1
if(isParallel){
    for(size_t i = 0;i < param.Den_thread_size;i++){
        sprintf( filename, "DenPartfig%d-%d-%d.plt", currentLevel, 0,i);
        ac_Den[i]->args_mynlp->m_pDB->OutputGnuplotFigureWithSpecifiedModules(filename,ac_Den[i]->args_mynlp->modules);
    }
    for(size_t i = 0;i < param.WL_thread_size;i++){
        sprintf( filename, "WLPartfig%d-%d-%d.plt", currentLevel, 0,i);
        ac_WL[i]->args_mynlp->m_pDB->OutputGnuplotFigureWithSpecifiedModules(filename,ac_WL[i]->args_mynlp->modules);
    }
}
#endif

    printf("\n======================================================================================\n");
    printf("cg begin\n");
    //lock      = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_init(&lock,nullptr);
    cond_lock = PTHREAD_COND_INITIALIZER;
    Xmerge = new double[size_module_2];
    //save_x_threads.resize(tot_tasks);
    for(int ite = 0;ite < maxIte;ite++){
     //   if( _weightWire <= 0.001) break;
        int innerIte = 0;
        double old_obj = DBL_MAX;
        double last_obj_value = DBL_MAX;

        memcpy(x_bak,x,sizeof(double)*size_module_2);

        while( true ){
            innerIte++;

            if( innerIte % checkStep == 0){
                printf(".");
                fflush(stdout);
                old_obj = last_obj_value;
                eval_f(size_module_2,x,_expX,true,obj_value);
                last_obj_value = obj_value;
            }
            if( innerIte % checkStep == 0 ){

                if( innerIte % (2*checkStep) == 0){
                    UpdateBlockPosition( x );
                    if( m_pDB->CalcHPWL() > bestLegalWL ){
                        printf("best legal WL\n");
                        fflush(stdout);
                        break;
                    }
                }

                UpdateDensityGrid(size_module_2,x);
                totalOverDen = GetTotalOverDensity();
                totalOverDenLB = GetTotalOverDensityLB();
                totalOverPotential = GetTotalOverPotential();

                lastTotalOver = over;
                over = min( totalOverPotential, totalOverDen ); // TEST

                if( !startDecreasing
                        && over < lastTotalOver
                        && ite >= 1
                        && innerIte >= 6)
                {
                    printf( ">>" );
                    fflush( stdout );
                    startDecreasing = true;
                }

                // 2005-03-11: meet the constraint
                if( startDecreasing && over < target_density ){
                    printf("over < target_density");
                    fflush(stdout);
                    break;
                }
                if( obj_value >= param.precision * old_obj){
                    printf("obj_value stop decreasing");
                    fflush(stdout);
                    break;
                }
            }
            //  update modules positions
            if(!isParallel){// no parallel
                if(!GD){
                    swap( last_grad_f, grad_f );
                    eval_grad_f( size_module_2, x, _expX, true, grad_f );
                    AdjustForce( size_module_2, x, grad_f );
                    if(innerIte == 1){
                        for(size_t i = 0;i < size_module_2;i++)
                            grad_f[i] = -grad_f[i];
                    }
                    else{
                        double beta;
                        FindBeta( size_module_2, grad_f, last_grad_f, beta );
                        for( size_t i=0; i<size_module_2; i++ )
                            grad_f[i] = -grad_f[i] + beta * last_grad_f[i];
                    }
                    // Calculate a_k (step size)
                   LineSearch( size_module_2, x, grad_f, stepSize );
                    // stepSize = param.step;
                    // Update X. (x_{k+1} = x_{k} + \alpha_k * d_k)
                    double move;
                    for( size_t i=0; i<size_module_2; i++ ){
                        move = grad_f[i] * stepSize;
                        x[i] += move;
                    }
                }
                else{// using GD
                    memcpy(xBak,x,sizeof(double)*size_module_2);
                    if(0){
                        eval_grad_f( size_module_2, x, _expX, true, grad_f );
                        DensityGrad(grad_density,whole_nets,bin);
                        WireLengthGrad(grad_wirelength,whole_nets);
                        double sum_diff = 0;
                        for(size_t i = 0;i < size_module_2;i++){
                            double temp = grad_density[i]*1 +
                                    grad_wirelength[i]*1;
                            double delta =  grad_f[i] - temp;
                            printf("gradx = %.4f, grad_potential = %.4f, grad_wire = %.4f\n",
                                   grad_f[i],grad_potential[i]*_weightDensity,grad_wire[i]*_weightWire);
                            printf("gradx = %.4f, grad_potential = %.4f, grad_wire = %.4f\ndelta = %.4f\n",
                                   temp,grad_density[i],
                                   grad_wirelength[i],
                                   delta);
                            sum_diff += delta;
                       //     getchar();
                        }
                        printf("sum_diff = %4f\n",sum_diff);
                        getchar();
                    }
                    else{


                        DensityGradByModule(grad_density,modules,bin);
                        WireLengthGradByModule(grad_wirelength,modules);

                        for(size_t i = 0;i < size_module_2;i++){
                            grad[i] = grad_density[i] + grad_wirelength[i];
                        }
                        //DensityAdjustForce(size_module_2,grad,whole_nets);
                        //AdjustForce(size_module_2,x,grad_f);
                        LineSearch(size_module_2,x,grad,stepsize);// 0.000687    0.005845
                        //stepsize = param.step;                       
                        for(size_t i = 0;i < size_module_2/2;i++){
                            x[2*i]   -= stepsize*grad[2*i];
                            x[2*i+1] -= stepsize*grad[2*i+1];
                        }
                    }


                }

            }
            else{// parallel
               // innerIte += (checkStep - 1);
                job_done = 0;
                done_flag = false;
#if 1
                DensityGradByModule(grad_density,modules,bin);
                WireLengthGradByModule(grad_wirelength,modules);
                for(size_t i = 0;i < size_module_2;i++){
                      grad[i] = grad_density[i] + grad_wirelength[i];
                }
                //AdjustForce(size_module_2,x,grad_f);
                LineSearch(size_module_2,x,grad,stepsize);
                stepsize *= parallelScale;
#endif
                vector<double> grad_WL_parallel;
                vector<double> grad_Density_parallel;
                //stepsize = param.step;
                //stepsize = (currentLevel == 1) ? 0.000006507 :  0.005;
                memcpy(xBak,x,sizeof(double)*size_module_2);
                for(unsigned int i_thread = 0;i_thread < param.Den_thread_size;i_thread++){
                     // for density
                    memcpy(ac_Den[i_thread]->args_mynlp->x,x,sizeof(double)*size_module_2);
                    ac_Den[i_thread]->stepsize = stepsize;//0.000967    0.009235
                    ac_Den[i_thread]->args_mynlp->_weightWire = _weightWire;
                    ac_Den[i_thread]->args_mynlp->_weightDensity = _weightDensity;
                    assert( threadpool_add(pool[i_thread],MyNLP::task_Density,(void*)(ac_Den[i_thread]),0) == 0 );
                }
                for(unsigned int i_thread = 0;i_thread < param.WL_thread_size;i_thread++){
                    // for wirelength
                    memcpy(ac_WL[i_thread]->args_mynlp->x,x,sizeof(double)*size_module_2);
                    ac_WL[i_thread]->stepsize = stepsize;
                    ac_WL[i_thread]->args_mynlp->_weightWire = _weightWire;
                    ac_WL[i_thread]->args_mynlp->_weightDensity = _weightDensity;
                    assert( threadpool_add(pool[i_thread + param.Den_thread_size],MyNLP::task_WireLength,(void*)(ac_WL[i_thread]),0) == 0 );
                }
                if(param.debug){
                    stepsize /= parallelScale;
                    for(size_t i = 0;i < size_module_2/2;i++){
                        x[2*i]   -= stepsize*grad[2*i];
                        x[2*i+1] -= stepsize*grad[2*i+1];
                    }
                    BoundX(size_module_2,x,x_l,x_u);
                }

                pthread_mutex_lock(&lock);
                pthread_cond_wait(&cond_lock,&lock);
                pthread_mutex_unlock(&lock);
                /*while(1){
                    pthread_mutex_lock(&lock);
                    if(job_done == tot_tasks){
                        pthread_mutex_unlock(&lock);
                        break;
                    }
                    pthread_mutex_unlock(&lock);
                }*/
            }


            double xdiff = 0;
            double *x_bak_bak;
            if(param.debug){
                x_bak_bak = new double[size_module_2];
                for(size_t i = 0; i < size_module_2;i++){
                    x_bak_bak[i] = x[i];
                }
            }
            // merge
#if 1
            if(isParallel){
                for(size_t i = 0;i < size_module_2;i++){
                    x[i] = Xmerge[i]/tot_tasks;
                    if(param.debug){
                        double delta = fabs(Xmerge[i] - (ac_Den[0]->args_mynlp->x[i] + ac_WL[0]->args_mynlp->x[i]));
                        assert(delta <= 0.00000001);
                        //printf("%4.3f\n",Xmerge[i] - (ac_Den[0]->args_mynlp->x[i] + ac_WL[0]->args_mynlp->x[i]));
                        xdiff += fabs(x[i] - x_bak_bak[i]);
                        //printf("%4.3f\n",fabs(x[i] - x_bak_bak[i]));
                    }

                    Xmerge[i] = 0;
                }
            }
            if(param.debug){
                //printf("\ntot diff = %4.3f,mean diff = %4.3f\n",xdiff,xdiff/size_module_2);
                param.xdiff.push_back(xdiff);
                //getchar();
            }
#endif
            if(0){
                double *x_old = new double[size_module_2];
                memcpy(x_old,x,sizeof(double)*size_module_2);
                double *expX_old = new double[size_module_2];
                for(size_t i = 0;i < size_module_2;i++){
                    assert(x[i] == x_old[i]);
                }
                BoundX( size_module_2, x_old, x_l, x_u );
                DensityBoundX(size_module_2);

                UpdateExpValueForEachCell( size_module_2, x, expX_old, _alpha );
                WireLengthUpdateExpValueForCell(whole_nets,size_module_2);
                double sum1 = 0,sum2 = 0;
                for(size_t i = 0;i < size_module_2;i++){
                    double delta1 = x_old[i] - x[i];
                    double delta2 = expX_old[i] - _expX[i];
                    printf("oldx[%5d] = %4.3f,newx[%5d] = %4.3f,delta = %4.3f\n",
                           i,x_old[i],i,x[i],delta1);
                //    printf("oldexpX[%5d] = %4.3f,newexpX[%5d] = %4.3f,delta = %4.3f\n",\
                           i,expX_old[i],i,_expX[i],delta2);
                    sum1 += fabs(delta1);
                    sum2 += fabs(delta2);
                }
                printf("sum1 = %4.3f,sum2 = %4.3f",sum1,sum2);
                getchar();
            }

            if(1){//

                BoundX( size_module_2, x, x_l, x_u );
                UpdateExpValueForEachCell( size_module_2, x, _expX, _alpha );
                UpdateExpValueForEachPin( size_module_2, x, _expPins, _alpha );
                UpdateNetsSumExp( x, _expX );
                UpdatePotentialGrid(x);
            }
            else{  //
                DensityBoundX(size_module_2);
                UpdatePotentialGrid(x);
                //DensityUpdatePotentialGridByModule(modules,bin);
                //WireLengthUpdateExpValueForCell(whole_nets,size_module_2);
                WireLengthUpdateExpValueForCellByModule(modules,size_module_2);
                WireLengthUpdateExpValueForPin(whole_nets);
                WireLengthUpdateNetsSumExp(whole_nets);
            }
#if 0
            if(innerIte == 4){
                nets_cluster_Den.clear();
                nets_cluster_Den.resize(bin_size);
                for(size_t i = 0;i < bin_size;i++)
                    nets_cluster_Den[i].resize(bin_size);

                for(unsigned int i = 0;i < m_pDB->m_nets.size();i++){
                    double c_x,c_y;
                    get_center_position(x,c_x,c_y,i,0);//0: ortho_center,1: mean_center
                    unsigned int n_x = (c_x - m_pDB->m_coreRgn.left)/(bucket_width );
                    unsigned int n_y = (c_y - m_pDB->m_coreRgn.bottom)/(bucket_height);
                    nets_cluster_Den[n_x][n_y].push_back(i);
                }
                for(size_t i = 0 ;i < bin_size;i++)
                    for(size_t j = 0;j < bin_size;j++){// re-arrange nets
                        size_t index = i*bin_size + j;
                        ac_Den[index]->args_snc_Id.clear();
                        ac_Den[index]->args_snc_Id.assign(nets_cluster_Den[i][j].begin(),nets_cluster_Den[i][j].end());
                        GetModulesByNets(ac_Den[index]->args_snc_Id,ac_Den[index]->args_mynlp->modules);

                        ac_WL[index]->args_snc_Id.clear();
                        ac_WL[index]->args_snc_Id.assign(nets_cluster_Den[i][j].begin(),nets_cluster_Den[i][j].end());
                        GetModulesByNets(ac_WL[index]->args_snc_Id,ac_WL[index]->args_mynlp->modules);
                    }
            }
#endif
            save_x_threads.clear();
            printf("*");
            fflush(stdout);
        } // end while
        if(param.debug){
        printf("\n");
            for(size_t i = 0;i < param.xdiff.size();i++)
                printf("%4.3f\t",param.xdiff[i]);
            param.xdiff.clear();
            getchar();
        }

        totalIte += innerIte;

        UpdateDensityGrid(size_module_2,x);
        UpdateBlockPosition(x);

        //        double nnb_real = GetNonZeroDensityGridPercent();
        maxDen = GetMaxDensity();
        totalOverDen = GetTotalOverDensity();
        totalOverDenLB = GetTotalOverDensityLB();
        totalOverPotential = GetTotalOverPotential();
        over = min( totalOverPotential, totalOverDen ); // TEST
        if( param.bShow )
        {
            if(ite+1 < 10){
                sprintf( filename, "PLfig%d-0%d.plt", currentLevel, ite+1 );
                m_pDB->OutputGnuplotFigure( filename, false, false );
            }
            else{
                sprintf( filename, "PLfig%d-%d.plt", currentLevel, ite+1 );
                m_pDB->OutputGnuplotFigure( filename, false, false );
            }
            if(isParallel){
                for(size_t i = 0;i < param.Den_thread_size;i++){
                    sprintf( filename, "DenPartfig%d-%d-%d.plt", currentLevel, ite+1,i);
                    ac_Den[i]->args_mynlp->m_pDB->OutputGnuplotFigureWithSpecifiedModules(filename,ac_Den[i]->args_mynlp->modules);
                }
                for(size_t i = 0;i < param.WL_thread_size;i++){
                    sprintf( filename, "WLPartfig%d-%d-%d.plt", currentLevel, ite+1,i);
                    ac_WL[i]->args_mynlp->m_pDB->OutputGnuplotFigureWithSpecifiedModules(filename,ac_WL[i]->args_mynlp->modules);
                }
            }
        }

        if(param.bShow){
            eval_f(size_module_2,x,_expX,true,obj_value);
            printf( "\n%d-%2d HPWL= %.0f  maxDen= %.2f OD= %.4f ODLB=%.4f OP=%.4f  LTime= %.1fm  objValue=%.0f, WireW= %.0f,totalWL = %.4f,DenW= %.4f,density = %.4f \n",
                    currentLevel, ite+1, m_pDB->CalcHPWL(),
                    maxDen, totalOverDen, totalOverDenLB, totalOverPotential,
                    double(seconds()-time_start)/60.0,
                    obj_value,
                    _weightWire,totalWL,
                    _weightDensity,density
                    );
        }

        if (1 && ite >= 2 && m_lookAheadLegalization && over < target_density + 0.10 )//2
        {
            //printf("in the loop,over < target_density + 0.1\n");
        //    UpdateBlockPosition(x);
            double hpwl = m_pDB->CalcHPWL();
            if (hpwl > bestLegalWL)
            {
                printf("Stop. Good enough 1.\n");
                break;
            }

            lookAheadLegalCount++;
            double oldWL = hpwl;
            CTetrisLegal legal(*m_pDB);

            double scale = 0.85;
            if (givenTargetUtil < 1.0 && givenTargetUtil > 0)
            {
                scale = 0.9;
            }

            double legalStart = seconds();
            bool   bLegal = legal.Solve( givenTargetUtil, false, false, scale);
            double legalTime = seconds() - legalStart;

            totalLegalTime += legalTime;

            if (param.bShow)
            {
                printf("LAL Time: %.2f\n", legalTime);
            }

            if (bLegal)
            {
                double WL = m_pDB->GetHPWLdensity(givenTargetUtil);
                if (param.bShow)
                {
                    m_pDB->ShowDensityInfo();
                }

                if (WL < bestLegalWL)
                {
                    LALnoGoodCount = 0;
                    if (param.bShow)
                    {
                        printf("SAVE BEST! (HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n",
                               m_pDB->GetHPWLp2p(), WL, (WL - oldWL) / oldWL * 100);
                    }
                    bestLegalWL = WL;
                    hasBestLegalSol = true;

                    for (size_t i = 0; i <  m_pDB->m_modules.size(); i++)
                    {
                        xBest[2 * i]     = m_pDB->m_modules[i].m_cx;
                        xBest[2 * i + 1] = m_pDB->m_modules[i].m_cy;
                        //                        x[2 * i]     = m_pDB->m_modules[i].m_cx;
                        //                        x[2 * i + 1] = m_pDB->m_modules[i].m_cy;

                    }
                }
                else {
                    if( param.bShow )
                        printf( "(HPWL=%.0f)(dHPWL= %.0f)(%.2f%%)\n", m_pDB->GetHPWLp2p(), WL, (WL-oldWL)/oldWL*100 );
                    if( (WL-oldWL)/oldWL < 0.075 ) // 0.01
                    {
                        if( param.bShow )
                            printf( "Stop. Good enough 2\n" );
                        break;
                    }
                    LALnoGoodCount++;
                    if( LALnoGoodCount >= 2 ){
                        printf("LALnoGoodCount >= 2\n");
                        break;
                    }
                }
            }
        } //end of m_lookAheadLegalization
        if (ite >= 2)//2  hdgao changed
        {
            if ( startDecreasing && over < target_density)// + 0.01
            {
                printf("Meet constraint!\n");
                break;
            }


            if (ite > 3 && totalOverPotential > lastTotalOverPotential &&
                    totalOverPotential < 1.4)
            {
                printf("Cannot further reduce!\n");
                break;
            }
        }


       // _weightWire /= m_incFactor;   //  _weightDensity *=2;
        _weightDensity *=2;
        lastTotalOverPotential = totalOverPotential;
#if 0
        if(ite == 3){
            nets_cluster_Den.clear();
            nets_cluster_Den.resize(bin_size);
            for(size_t i = 0;i < bin_size;i++)
                nets_cluster_Den[i].resize(bin_size);

            for(unsigned int i = 0;i < m_pDB->m_nets.size();i++){
                double c_x,c_y;
                get_center_position(x,c_x,c_y,i,0);//0: ortho_center,1: mean_center
                unsigned int n_x = (c_x - m_pDB->m_coreRgn.left)/(bucket_width );
                unsigned int n_y = (c_y - m_pDB->m_coreRgn.bottom)/(bucket_height);
                nets_cluster_Den[n_x][n_y].push_back(i);
            }
            for(size_t i = 0 ;i < bin_size;i++)
                for(size_t j = 0;j < bin_size;j++){// re-arrange nets
                    size_t index = i*bin_size + j;
                    ac_Den[index]->args_snc_Id.clear();
                    ac_Den[index]->args_snc_Id.assign(nets_cluster_Den[i][j].begin(),nets_cluster_Den[i][j].end());
                    GetModulesByNets(ac_Den[index]->args_snc_Id,ac_Den[index]->args_mynlp->modules);

                    ac_WL[index]->args_snc_Id.clear();
                    ac_WL[index]->args_snc_Id.assign(nets_cluster_Den[i][j].begin(),nets_cluster_Den[i][j].end());
                    GetModulesByNets(ac_WL[index]->args_snc_Id,ac_WL[index]->args_mynlp->modules);
                }
        }
#endif
    }// end of "for (ite = 0; ite < maxIte; ite++)"

    if (hasBestLegalSol)
    {
        memcpy(x, xBest, sizeof(double) * size_module_2);
    }

    UpdateBlockPosition(x);

    if( lookAheadLegalCount > 0 && param.bShow )
    {
        printf( "LAL: Total Count: %d\n", lookAheadLegalCount );
        printf( "LAL: Total CPU: %.2f\n", totalLegalTime );
        sprintf( filename, "util-global.dat" );
        CPlaceBin placeBin( *m_pDB );
        placeBin.CreateGrid( m_pDB->m_rowHeight * 10.0 );
        placeBin.OutputBinUtil( filename );
    }
    printf("=======================avg stepsize = %f\n",param.sum_step/param.count);
    param.stepRecord.clear();
    param.sum_step = 0;
    param.count = 0;
    delete [] grad_f;
    delete [] last_grad_f;	// for computing conjugate gradient direction
    return hasBestLegalSol;
}

/*
void MyNLP::mynlp_copy(MyNLP* source){
    // just copy what we want
    target->_cellPotentialNorm.resize( m_pDB->m_modules.size() );

    target->x        = new double [ 2 * m_pDB->m_modules.size() ];
    target->_expX    = new double [ 2 * m_pDB->m_modules.size() ];
    target->_expPins = new double [ 2 * m_pDB->m_pins.size() ];
    target->x_l      = new double [ 2 * m_pDB->m_modules.size() ];
    target->x_u      = new double [ 2 * m_pDB->m_modules.size() ];

    target->m_usePin.resize( m_pDB->m_modules.size() );
    target->SetUsePin();

    target->m_nets_sum_exp_xi_over_alpha.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_exp_yi_over_alpha.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_exp_inv_xi_over_alpha.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_exp_inv_yi_over_alpha.resize( m_pDB->m_nets.size(), 0 );

    target->m_nets_sum_p_x_pos.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_p_y_pos.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_p_inv_x_pos.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_p_inv_y_pos.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_p_x_neg.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_p_y_neg.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_p_inv_x_neg.resize( m_pDB->m_nets.size(), 0 );
    target->m_nets_sum_p_inv_y_neg.resize( m_pDB->m_nets.size(), 0 );

    target->grad_wire.resize( 2 * m_pDB->m_modules.size(), 0.0 );
    target->grad_potential.resize( 2 * m_pDB->m_modules.size(), 0.0 );

    //flpengchange 2015 12.04
    target->sub_grad_wire.resize(2 * m_pDB->m_modules.size(), 0.0);
    target->sub_grad_potential.resize(2 * m_pDB->m_modules.size(), 0.0);


    // hdgao added
    target->part_grad_wire.resize(2 * m_pDB->m_modules.size(), 0.0);
    target->part_grad_potential.resize(2 * m_pDB->m_modules.size(), 0.0);
    target->part_grad.resize(2*m_pDB->m_modules.size(),0.0);
}
*/
void MyNLP::CG(MyNLP* cur_mynlp,int size_x, double* x_bak, bool is_2D){
    double step_sum = 0;
    int n_step_sum = 0;
    unsigned int size_nets_cluster = cur_mynlp->nets_cluster.size();
    for(unsigned int i_nets_cluster = 0; i_nets_cluster < size_nets_cluster ;i_nets_cluster++){

        vector< double > last_grad_f,grad_f;
        last_grad_f.resize(size_x,0);
        grad_f.resize(size_x,0);
        vector<int> RegionModule;
        RegionModule.clear();
        //              RegionModule.assign(nets_black_modules[i_nets_cluster].begin(),nets_black_modules[i_nets_cluster].end());
        int cg_numIte = 2;
        for(int i_cg = 0;i_cg < cg_numIte;i_cg++){// 5 -> 3
            last_grad_f.assign( grad_f.begin(),grad_f.end() );
            if(i_cg == 0)  cur_mynlp->updateBlackBounds(i_nets_cluster);
            cur_mynlp->smooth_part_eval_grad_f(cur_mynlp->x,cur_mynlp->nets_cluster_modules_bounds[i_nets_cluster],
                                               cur_mynlp->nets_black_modules[i_nets_cluster],
                                               RegionModule,grad_f); // need change


            if(is_2D)
                cur_mynlp->PartAdjustForce_2D(grad_f,RegionModule);
            else
                cur_mynlp->PartAdjustForce(grad_f,RegionModule);
            //calculate d_k
            double beta1,beta2,beta;
            if(i_cg == 0){// first
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    size_t moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1];
                }
            }
            else{
                if(is_2D)
                    cur_mynlp->PartFindBeta_2D(grad_f,last_grad_f,RegionModule,beta1,beta2);
                else{
                    cur_mynlp->PartFindBeta(grad_f,last_grad_f,RegionModule,beta);
                    beta1 = beta2 = beta;
                }
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ]    + beta1 * last_grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1] + beta2 * last_grad_f[ 2*moduleId + 1];
                }
            }

            //calculate a_k , the step size
            double stepSize1,stepSize2;
            if(is_2D){
                cur_mynlp->PartLineSearch_2D(grad_f,RegionModule,stepSize1,stepSize2);
                step_sum = step_sum + ( stepSize1+ stepSize2)/2;
                n_step_sum ++;
            }
            else{
                cur_mynlp->PartLineSearch(grad_f,RegionModule,stepSize);
                stepSize1 = stepSize2 = stepSize;
                step_sum += stepSize;
                n_step_sum ++;
            }
            for (size_t i = 0; i < RegionModule.size(); i++){
                int moduleId = RegionModule[i];
                cur_mynlp->x[2*moduleId]   += grad_f[2*moduleId]  * stepSize1;
                cur_mynlp->x[2*moduleId+1] += grad_f[2*moduleId+1]* stepSize2;
            }

            cur_mynlp->BoundX(size_x, cur_mynlp->x, cur_mynlp->x_l, cur_mynlp->x_u);
            double time_used = seconds();
            vector<int> right_mId;
            for(int i = 0; i < size_x/2;i++){
                if( x_bak[2*i] != cur_mynlp->x[2*i] || x_bak[2*i + 1] != cur_mynlp->x[2*i + 1])
                    right_mId.push_back(i);
            }
            cur_mynlp->part_UpdatePotentialGrid(cur_mynlp->x,x_bak,right_mId,right_mId.size());// to update m_potential

            for(int i = 0; i < right_mId.size();i++){
                x_bak[ 2*right_mId[i] ]   = cur_mynlp->x[ 2*right_mId[i] ];
                x_bak[ 2*right_mId[i]+1 ] = cur_mynlp->x[ 2*right_mId[i]+1 ];
            }

            // time-consuming here

            cur_mynlp->UpdateExpValueForEachCell( size_x, cur_mynlp->x, cur_mynlp->_expX, cur_mynlp->_alpha );  // to update _expx
            cur_mynlp->UpdateExpValueForEachPin( size_x, cur_mynlp->x, cur_mynlp->_expPins, cur_mynlp->_alpha );// to update _expPins
            cur_mynlp->UpdateNetsSumExp( cur_mynlp->x, cur_mynlp->_expX );// to update
            time_grad_wl += seconds() - time_used;
            //             UpdatePotentialGrid(x);
        }// end cg
    }// end traversal nets
    cur_mynlp->PL_step = step_sum/n_step_sum;
    //  printf("================================%f===============================\n",step_sum/n_step_sum);
}
void MyNLP::CG(MyNLP* cur_mynlp,int size_x, double *x_bak, bool is_2D,// process some specified nets cluster of nets_cluster
               const vector<int>& specified_nets_cluster){
    unsigned int size_nets_cluster = specified_nets_cluster.size();
    double *x_bakbb = new double(size_x);
    memcpy(x_bakbb,cur_mynlp->x,sizeof(double)*size_x);
    for(unsigned int i_nets_cluster = 0; i_nets_cluster < size_nets_cluster ;i_nets_cluster++){

        int snc = specified_nets_cluster[i_nets_cluster];
        vector< double > last_grad_f,grad_f;
        last_grad_f.resize(size_x,0);
        grad_f.resize(size_x,0);
        vector<int> RegionModule;
        RegionModule.clear();
        //              RegionModule.assign(nets_black_modules[i_nets_cluster].begin(),nets_black_modules[i_nets_cluster].end());
        int cg_numIte = 5;
        for(int i_cg = 0;i_cg < cg_numIte;i_cg++){// 5 -> 3
            last_grad_f.assign( grad_f.begin(),grad_f.end() );
            if(i_cg == 0)  cur_mynlp->updateBlackBounds(snc);
            cur_mynlp->smooth_part_eval_grad_f(cur_mynlp->x,cur_mynlp->nets_cluster_modules_bounds[snc],
                                               cur_mynlp->nets_black_modules[snc],
                                               RegionModule,grad_f); // need change


            if(is_2D)
                cur_mynlp->PartAdjustForce_2D(grad_f,RegionModule);
            else
                cur_mynlp->PartAdjustForce(grad_f,RegionModule);
            //calculate d_k
            double beta1,beta2,beta;
            if(i_cg == 0){// first
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1];
                }
            }
            else{
                if(is_2D)
                    cur_mynlp->PartFindBeta_2D(grad_f,last_grad_f,RegionModule,beta1,beta2);
                else{
                    cur_mynlp->PartFindBeta(grad_f,last_grad_f,RegionModule,beta);
                    beta1 = beta2 = beta;
                }
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ]    + beta1 * last_grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1] + beta2 * last_grad_f[ 2*moduleId + 1];
                }
            }

            //calculate a_k , the step size
            double stepSize1,stepSize2;
            if(is_2D)
                cur_mynlp->PartLineSearch_2D(grad_f,RegionModule,stepSize1,stepSize2);
            else{
                cur_mynlp->PartLineSearch(grad_f,RegionModule,stepSize);
                stepSize1 = stepSize2 = stepSize;
            }
            for (size_t i = 0; i < RegionModule.size(); i++){
                int moduleId = RegionModule[i];
                cur_mynlp->x[2*moduleId]   += grad_f[2*moduleId]  * stepSize1;
                cur_mynlp->x[2*moduleId+1] += grad_f[2*moduleId+1]* stepSize2;
            }

            cur_mynlp->BoundX(size_x, cur_mynlp->x, cur_mynlp->x_l, cur_mynlp->x_u);
            double time_used = seconds();
            vector<int> right_mId;
            for(int i = 0; i < size_x/2;i++){
                if( x_bakbb[2*i] != cur_mynlp->x[2*i] || x_bakbb[2*i + 1] != cur_mynlp->x[2*i + 1])
                    right_mId.push_back(i);
            }
            cur_mynlp->part_UpdatePotentialGrid(cur_mynlp->x,x_bakbb,right_mId,right_mId.size());// to update m_potential

            for(int i = 0; i < right_mId.size();i++){
                x_bakbb[ 2*right_mId[i] ]   = cur_mynlp->x[ 2*right_mId[i] ];
                x_bakbb[ 2*right_mId[i]+1 ] = cur_mynlp->x[ 2*right_mId[i]+1 ];
            }

            // time-consuming here

            cur_mynlp->UpdateExpValueForEachCell( size_x, cur_mynlp->x, cur_mynlp->_expX, cur_mynlp->_alpha );  // to update _expx
            cur_mynlp->UpdateExpValueForEachPin( size_x, cur_mynlp->x, cur_mynlp->_expPins, cur_mynlp->_alpha );// to update _expPins
            cur_mynlp->UpdateNetsSumExp( cur_mynlp->x, cur_mynlp->_expX );// to update
            time_grad_wl += seconds() - time_used;
            //             UpdatePotentialGrid(x);
        }// end cg
    }// end traversal nets
}


void MyNLP::CG(MyNLP* cur_mynlp, int size_x, double *x_bak, bool is_2D, // process some specified nets cluster of nets_cluster
               const vector< vector<int> >& cur_nets_cluster,
               const vector< vector<int> >& cur_nets_black_modules,
               vector<vector<int> > cur_nets_cluster_modules_bounds,
               const vector<int>& specified_nets_cluster){
    double step_sums = 0;
    int n_step_sums = 0;
    unsigned int size_nets_cluster = specified_nets_cluster.size();
    for(unsigned int i_nets_cluster = 0; i_nets_cluster < size_nets_cluster ;i_nets_cluster++){

        int snc = specified_nets_cluster[i_nets_cluster];
        vector< double > last_grad_f,grad_f;
        last_grad_f.resize(size_x,0);
        grad_f.resize(size_x,0);
        vector<int> RegionModule;
        RegionModule.clear();
        //              RegionModule.assign(nets_black_modules[i_nets_cluster].begin(),nets_black_modules[i_nets_cluster].end());
        int cg_numIte = 5;
        for(int i_cg = 0;i_cg < cg_numIte;i_cg++){// 5 -> 3
            last_grad_f.assign( grad_f.begin(),grad_f.end() );
            if(i_cg == 0)  cur_mynlp->updateBlackBounds(snc,cur_nets_cluster,cur_nets_black_modules,cur_nets_cluster_modules_bounds);
            cur_mynlp->smooth_part_eval_grad_f(cur_mynlp->x,cur_nets_cluster_modules_bounds[snc],
                                               cur_nets_black_modules[snc],
                                               RegionModule,grad_f); // need change


            if(is_2D)
                cur_mynlp->PartAdjustForce_2D(grad_f,RegionModule);
            else
                cur_mynlp->PartAdjustForce(grad_f,RegionModule);
            //calculate d_k
            double beta1,beta2,beta;
            if(i_cg == 0){// first
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1];
                }
            }
            else{
                if(is_2D)
                    cur_mynlp->PartFindBeta_2D(grad_f,last_grad_f,RegionModule,beta1,beta2);
                else{
                    cur_mynlp->PartFindBeta(grad_f,last_grad_f,RegionModule,beta);
                    beta1 = beta2 = beta;
                }
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ]    + beta1 * last_grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1] + beta2 * last_grad_f[ 2*moduleId + 1];
                }
            }

            //calculate a_k , the step size
            double stepSize1,stepSize2;
            if(is_2D){
                cur_mynlp->PartLineSearch_2D(grad_f,RegionModule,stepSize1,stepSize2);
                step_sums = step_sums + (stepSize1 + stepSize2)/2;
                n_step_sums ++;
            }
            else{
                cur_mynlp->PartLineSearch(grad_f,RegionModule,stepSize);
                stepSize1 = stepSize2 = stepSize;
                step_sums += stepSize;
                n_step_sums ++;
            }
            // fixed step
            if(parallel_flag) stepSize1 = stepSize2 = stepSize = cur_mynlp->PL_step;
            for (size_t i = 0; i < RegionModule.size(); i++){
                int moduleId = RegionModule[i];
                cur_mynlp->x[2*moduleId]   += grad_f[2*moduleId]  * stepSize1;
                cur_mynlp->x[2*moduleId+1] += grad_f[2*moduleId+1]* stepSize2;
            }

            cur_mynlp->BoundX(size_x, cur_mynlp->x, cur_mynlp->x_l, cur_mynlp->x_u);
            double time_used = seconds();
            vector<int> right_mId;
            for(int i = 0; i < size_x/2;i++){
                if( x_bak[2*i] != cur_mynlp->x[2*i] || x_bak[2*i + 1] != cur_mynlp->x[2*i + 1])
                    right_mId.push_back(i);
            }
            cur_mynlp->part_UpdatePotentialGrid(cur_mynlp->x,x_bak,right_mId,right_mId.size());// to update m_potential

            for(int i = 0; i < right_mId.size();i++){
                x_bak[ 2*right_mId[i] ]   = cur_mynlp->x[ 2*right_mId[i] ];
                x_bak[ 2*right_mId[i]+1 ] = cur_mynlp->x[ 2*right_mId[i]+1 ];
            }

            // time-consuming here

            cur_mynlp->UpdateExpValueForEachCell( size_x, cur_mynlp->x, cur_mynlp->_expX, cur_mynlp->_alpha );  // to update _expx
            cur_mynlp->UpdateExpValueForEachPin( size_x, cur_mynlp->x, cur_mynlp->_expPins, cur_mynlp->_alpha );// to update _expPins
            cur_mynlp->UpdateNetsSumExp( cur_mynlp->x, cur_mynlp->_expX );// to update
            time_grad_wl += seconds() - time_used;
            //             UpdatePotentialGrid(x);
        }// end cg
    }// end traversal nets
    //  printf("%f\n",step_sums/n_step_sums);
    //getchar();
}
void MyNLP::GD(MyNLP* cur_mynlp, int size_x, double *x_bak, bool is_2D, // process some specified nets cluster of nets_cluster
               const vector< vector<int> >& cur_nets_cluster,
               const vector< vector<int> >& cur_nets_black_modules,
               vector< vector<int> > cur_nets_cluster_modules_bounds,
               const vector<int>& specified_nets_cluster){
    for(int loop = 0; loop < 10 ; loop++){
        double stepSize,stepSize1,stepSize2;


        unsigned int size_nets_cluster = specified_nets_cluster.size();
        for(unsigned int i_nets_cluster = 0; i_nets_cluster < size_nets_cluster ;i_nets_cluster++){

            int snc = specified_nets_cluster[i_nets_cluster];
            vector< double > grad_f;
            grad_f.resize(size_x,0);
            vector<int> RegionModule;
            RegionModule.clear();
            //              RegionModule.assign(nets_black_modules[i_nets_cluster].begin(),nets_black_modules[i_nets_cluster].end());

            cur_mynlp->updateBlackBounds(snc,cur_nets_cluster,cur_nets_black_modules,cur_nets_cluster_modules_bounds);
            cur_mynlp->smooth_part_eval_grad_f(cur_mynlp->x,cur_nets_cluster_modules_bounds[snc],
                                               cur_nets_black_modules[snc],
                                               RegionModule,grad_f); // need change


            if(is_2D)
                cur_mynlp->PartAdjustForce_2D(grad_f,RegionModule);
            else
                cur_mynlp->PartAdjustForce(grad_f,RegionModule);

            //calculate a_k , the step size
            // fixed step
            if(parallel_flag) stepSize1 = stepSize2 = stepSize = 0.002;//cur_mynlp->PL_step;
            for (size_t i = 0; i < RegionModule.size(); i++){
                int moduleId = RegionModule[i];
                cur_mynlp->x[2*moduleId]   -= grad_f[2*moduleId]  * stepSize1;
                cur_mynlp->x[2*moduleId+1] -= grad_f[2*moduleId+1]* stepSize2;
            }

            cur_mynlp->BoundX(size_x, cur_mynlp->x, cur_mynlp->x_l, cur_mynlp->x_u);
            double time_used = seconds();
            vector<int> right_mId;
            for(int i = 0; i < size_x/2;i++){
                if( x_bak[2*i] != cur_mynlp->x[2*i] || x_bak[2*i + 1] != cur_mynlp->x[2*i + 1])
                    right_mId.push_back(i);
            }
            cur_mynlp->part_UpdatePotentialGrid(cur_mynlp->x,x_bak,right_mId,right_mId.size());// to update m_potential

            for(int i = 0; i < right_mId.size();i++){//update x_bak
                x_bak[ 2*right_mId[i] ]   = cur_mynlp->x[ 2*right_mId[i] ];
                x_bak[ 2*right_mId[i]+1 ] = cur_mynlp->x[ 2*right_mId[i]+1 ];
            }

            // time-consuming heregn

            cur_mynlp->UpdateExpValueForEachCell( size_x, cur_mynlp->x, cur_mynlp->_expX, cur_mynlp->_alpha );  // to update _expx
            cur_mynlp->UpdateExpValueForEachPin( size_x, cur_mynlp->x, cur_mynlp->_expPins, cur_mynlp->_alpha );// to update _expPins
            cur_mynlp->UpdateNetsSumExp( cur_mynlp->x, cur_mynlp->_expX );// to update
            time_grad_wl += seconds() - time_used;


            //             UpdatePotentialGrid(x);
        }// end traversal nets
        //   printf("%f\n",step_sums/n_step_sums);
        // getchar();
    }// end loop
}
void MyNLP::SGD_Density(MyNLP* cur_mynlp,size_t size_x,
                        const vector<int>& specified_nets_cluster,
                        const vector<int>& bin,
                        double stepsize,double maxloop){
 //   if(cur_mynlp->xBak == nullptr)
   //     cur_mynlp->xBak = new double[size_x];

    vector<double> grad;
    grad.resize(size_x,0);
    for(size_t loop = 0; loop < maxloop; loop++){

        //memcpy(cur_mynlp->xBak,cur_mynlp->x,sizeof(double)*size_x);
        for(size_t i = 0;i < size_x;i++)
            cur_mynlp->xBak[i] = cur_mynlp->x[i];

        //cur_mynlp->DensityGrad(grad,specified_nets_cluster,bin);
        cur_mynlp->DensityGradByModule(grad,cur_mynlp->modules,bin);
        //cur_mynlp->DensityAdjustForce(size_x,grad,specified_nets_cluster);
        for(size_t i = 0; i < size_x/2; i++){
            cur_mynlp->x[2*i]   -= stepsize*grad[2*i];
            cur_mynlp->x[2*i+1] -= stepsize*grad[2*i+1];
        }

        cur_mynlp->DensityBoundX(size_x);
        //cur_mynlp->DensityUpdatePotentialGrid(specified_nets_cluster,bin);
        cur_mynlp->DensityUpdatePotentialGridByModule(cur_mynlp->modules,bin);
        //cur_mynlp->UpdatePotentialGrid(cur_mynlp->x);
 //       cur_mynlp->DensityUpdateExpValueForCell();
 //       cur_mynlp->DensityUpdateExpValueForPin();
  //      cur_mynlp->DensityUpdateNetsSumExp();
    }
    //delete cur_mynlp->xBak;
}
void MyNLP::SGD_WireLength(MyNLP* cur_mynlp,size_t size_x,const vector<int>& specified_nets_cluster,
                        double stepsize,double maxloop){
    double SS = stepsize;
    vector<double> grad;
    grad.resize(size_x,0);
    for(size_t loop = 0; loop < maxloop; loop++){


   //     cur_mynlp->WireLengthGrad(grad,specified_nets_cluster);
        cur_mynlp->WireLengthGradByModule(grad,cur_mynlp->modules);
        //cur_mynlp->WireLengthAdjustForce(size_x,grad,specified_nets_cluster);
        for(size_t i = 0; i < size_x/2; i++){
            cur_mynlp->x[2*i]   -= ( grad[2*i]*SS );
            cur_mynlp->x[2*i+1] -= ( grad[2*i+1]*SS);
        }

        cur_mynlp->WireLengthBoundX(size_x);
//        cur_mynlp->WireLengthUpdatePotentialGrid();
//        cur_mynlp->WireLengthUpdateExpValueForCell(specified_nets_cluster,size_x);
        cur_mynlp->WireLengthUpdateExpValueForCellByModule(cur_mynlp->modules,size_x);
        cur_mynlp->WireLengthUpdateExpValueForPin(specified_nets_cluster);
        cur_mynlp->WireLengthUpdateNetsSumExp(specified_nets_cluster);
    }
}


void MyNLP::CG_mynlp(int size_x, double* x_bak, bool is_2D){
    unsigned int size_nets_cluster =nets_cluster.size();
    for(unsigned int i_nets_cluster = 0; i_nets_cluster < size_nets_cluster ;i_nets_cluster++){
        vector< double > last_grad_f,grad_f;
        last_grad_f.assign(size_x,0);
        grad_f.assign(size_x,0);
        vector<int> RegionModule;
        RegionModule.clear();
        //              RegionModule.assign(nets_black_modules[i_nets_cluster].begin(),nets_black_modules[i_nets_cluster].end());
        int cg_numIte = 5;
        for(int i_cg = 0;i_cg < cg_numIte;i_cg++){// 5 -> 3
            last_grad_f.assign( grad_f.begin(),grad_f.end() );
            if(i_cg == 0)  updateBlackBounds(i_nets_cluster);
            smooth_part_eval_grad_f(x,nets_cluster_modules_bounds[i_nets_cluster],
                                    nets_black_modules[i_nets_cluster],
                                    RegionModule,grad_f); // need change


            if(is_2D)
                PartAdjustForce_2D(grad_f,RegionModule);
            else
                PartAdjustForce(grad_f,RegionModule);
            //calculate d_k
            double beta1,beta2,beta;
            if(i_cg == 0){// first
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1];
                }
            }
            else{
                if(is_2D)
                    PartFindBeta_2D(grad_f,last_grad_f,RegionModule,beta1,beta2);
                else{
                    PartFindBeta(grad_f,last_grad_f,RegionModule,beta);
                    beta1 = beta2 = beta;
                }
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ]    + beta1 * last_grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1] + beta2 * last_grad_f[ 2*moduleId + 1];
                }
            }

            //calculate a_k , the step size
            double stepSize1,stepSize2;
            if(is_2D)
                PartLineSearch_2D(grad_f,RegionModule,stepSize1,stepSize2);
            else{
                PartLineSearch(grad_f,RegionModule,stepSize);
                stepSize1 = stepSize2 = stepSize;
            }
            for (size_t i = 0; i < RegionModule.size(); i++){
                int moduleId = RegionModule[i];
                x[2*moduleId]   += grad_f[2*moduleId]  * stepSize1;
                x[2*moduleId+1] += grad_f[2*moduleId+1]* stepSize2;
            }

            BoundX(size_x, x, x_l, x_u);
            double time_used = seconds();
            vector<int> right_mId;
            for(int i = 0; i < size_x/2;i++){
                if( x_bak[2*i] != x[2*i] || x_bak[2*i + 1] != x[2*i + 1])
                    right_mId.push_back(i);
            }
            part_UpdatePotentialGrid(x,x_bak,right_mId,right_mId.size());// to update m_potential

            for(int i = 0; i < right_mId.size();i++){
                x_bak[ 2*right_mId[i] ]   = x[ 2*right_mId[i] ];
                x_bak[ 2*right_mId[i]+1 ] = x[ 2*right_mId[i]+1 ];
            }

            // time-consuming here

            UpdateExpValueForEachCell( size_x, x, _expX, _alpha );  // to update _expx
            UpdateExpValueForEachPin( size_x, x, _expPins, _alpha );// to update _expPins
            UpdateNetsSumExp( x, _expX );// to update
            time_grad_wl += seconds() - time_used;
            //             UpdatePotentialGrid(x);
        }// end cg
    }// end traversal nets
}

void MyNLP::CG_mynlp(int size_x, double *x_bak, bool is_2D,// process some specified nets cluster of nets_cluster
                     //               const vector< vector<int> >& cur_nets_cluster,
                     //               const vector< vector<int> >& cur_nets_black_modules,
                     //               const vector< vector<int> >& cur_nets_cluster_modules_bounds,
                     const vector<int>& specified_nets_cluster){
    unsigned int size_nets_cluster = specified_nets_cluster.size();
    for(unsigned int i_nets_cluster = 0; i_nets_cluster < size_nets_cluster ;i_nets_cluster++){

        int snc = specified_nets_cluster[i_nets_cluster];
        vector< double > last_grad_f,grad_f;
        last_grad_f.resize(size_x,0);
        grad_f.resize(size_x,0);
        vector<int> RegionModule;
        RegionModule.clear();
        //              RegionModule.assign(nets_black_modules[i_nets_cluster].begin(),nets_black_modules[i_nets_cluster].end());
        int cg_numIte = 5;
        for(int i_cg = 0;i_cg < cg_numIte;i_cg++){// 5 -> 3
            last_grad_f.assign( grad_f.begin(),grad_f.end() );
            if(i_cg == 0)  updateBlackBounds(snc);
            smooth_part_eval_grad_f(x,nets_cluster_modules_bounds[snc],
                                    nets_black_modules[snc],
                                    RegionModule,grad_f); // need change


            if(is_2D)
                PartAdjustForce_2D(grad_f,RegionModule);
            else
                PartAdjustForce(grad_f,RegionModule);
            //calculate d_k
            double beta1,beta2,beta;
            if(i_cg == 0){// first
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1];
                }
            }
            else{
                if(is_2D)
                    PartFindBeta_2D(grad_f,last_grad_f,RegionModule,beta1,beta2);
                else{
                    PartFindBeta(grad_f,last_grad_f,RegionModule,beta);
                    beta1 = beta2 = beta;
                }
                for (unsigned int i = 0; i < RegionModule.size(); i++){
                    int moduleId = RegionModule[i];
                    grad_f[ 2*moduleId ]    = -grad_f[ 2*moduleId ]    + beta1 * last_grad_f[ 2*moduleId ];
                    grad_f[ 2*moduleId + 1] = -grad_f[ 2*moduleId + 1] + beta2 * last_grad_f[ 2*moduleId + 1];
                }
            }

            //calculate a_k , the step size
            double stepSize1,stepSize2;
            if(is_2D)
                PartLineSearch_2D(grad_f,RegionModule,stepSize1,stepSize2);
            else{
                PartLineSearch(grad_f,RegionModule,stepSize);
                stepSize1 = stepSize2 = stepSize;
            }
            for (size_t i = 0; i < RegionModule.size(); i++){
                int moduleId = RegionModule[i];
                x[2*moduleId]   += grad_f[2*moduleId]  * stepSize1;
                x[2*moduleId+1] += grad_f[2*moduleId+1]* stepSize2;
            }

            BoundX(size_x, x, x_l, x_u);
            double time_used = seconds();
            vector<int> right_mId;
            for(int i = 0; i < size_x/2;i++){
                if( x_bak[2*i] != x[2*i] || x_bak[2*i + 1] != x[2*i + 1])
                    right_mId.push_back(i);
            }
            part_UpdatePotentialGrid(x,x_bak,right_mId,right_mId.size());// to update m_potential

            for(int i = 0; i < right_mId.size();i++){
                x_bak[ 2*right_mId[i] ]   = x[ 2*right_mId[i] ];
                x_bak[ 2*right_mId[i]+1 ] = x[ 2*right_mId[i]+1 ];
            }

            // time-consuming here

            UpdateExpValueForEachCell( size_x, x, _expX, _alpha );  // to update _expx
            UpdateExpValueForEachPin( size_x, x, _expPins, _alpha );// to update _expPins
            UpdateNetsSumExp(x,_expX );// to update
            time_grad_wl += seconds() - time_used;
            //             UpdatePotentialGrid(x);
        }// end cg
    }// end traversal nets
}


void MyNLP::task_CG(void *args){
    args_class* AC = (args_class*)args;
    assert(AC != NULL);
    // update with x
    AC->args_mynlp->BoundX(AC->size_x, AC->args_mynlp->x,
                           AC->args_mynlp->x_l, AC->args_mynlp->x_u);
    AC->args_mynlp->UpdateExpValueForEachCell(AC->size_x, AC->args_mynlp->x,
                                              AC->args_mynlp->_expX, AC->args_mynlp->_alpha);
    AC->args_mynlp->UpdateExpValueForEachPin(AC->size_x, AC->args_mynlp->x,
                                             AC->args_mynlp->_expX, AC->args_mynlp->_alpha);
    AC->args_mynlp->UpdateNetsSumExp(AC->args_mynlp->x,
                                     AC->args_mynlp->_expX);
    AC->args_mynlp->UpdateBlockPosition(AC->args_mynlp->x);
    AC->args_mynlp->UpdateDensityGrid(AC->size_x,AC->args_mynlp->x);
    AC->args_mynlp->UpdatePotentialGrid(AC->args_mynlp->x);


    CG(AC->args_mynlp,AC->size_x,AC->x_bak,AC->is_2D,
       AC->args_mynlp->nets_cluster,
       AC->args_mynlp->nets_black_modules,
       AC->args_mynlp->nets_cluster_modules_bounds,
       AC->args_snc_Id);// need x back

    vector<double> temp;
    for(size_t i = 0;i < 2*AC->args_mynlp->m_pDB->m_nModules;i++)
        temp.push_back(AC->args_mynlp->x[i]);
    pthread_mutex_lock(&lock);
    save_x_threads.push_back(temp);
    job_done++;
    //  printf("thread %d done!\n",AC->thread_number);
    if(job_done == THREAD_size)
        done_flag = true;
    pthread_mutex_unlock(&lock);
}
void MyNLP::task_GD(void *args){
    args_class* AC = (args_class*)args;
    assert(AC != nullptr);

    //  printf("this is thread %d",AC->thread_number);
    // update with x
    AC->args_mynlp->BoundX(AC->size_x, AC->args_mynlp->x,
                           AC->args_mynlp->x_l, AC->args_mynlp->x_u);
    AC->args_mynlp->UpdateExpValueForEachCell(AC->size_x, AC->args_mynlp->x,
                                              AC->args_mynlp->_expX, AC->args_mynlp->_alpha);
    AC->args_mynlp->UpdateExpValueForEachPin(AC->size_x, AC->args_mynlp->x,
                                             AC->args_mynlp->_expX, AC->args_mynlp->_alpha);
    AC->args_mynlp->UpdateNetsSumExp(AC->args_mynlp->x,
                                     AC->args_mynlp->_expX);
    AC->args_mynlp->UpdateBlockPosition(AC->args_mynlp->x);
    AC->args_mynlp->UpdateDensityGrid(AC->size_x,AC->args_mynlp->x);
    AC->args_mynlp->UpdatePotentialGrid(AC->args_mynlp->x);


    GD(AC->args_mynlp,AC->size_x,AC->x_bak,AC->is_2D,
       AC->args_mynlp->nets_cluster,
       AC->args_mynlp->nets_black_modules,
       AC->args_mynlp->nets_cluster_modules_bounds,
       AC->args_snc_Id);// need x back

    vector<double> temp;
    for(size_t i = 0;i < 2*AC->args_mynlp->m_pDB->m_modules.size();i++)
        temp.push_back(AC->args_mynlp->x[i]);
    pthread_mutex_lock(&lock);
    save_x_threads.push_back(temp);
    job_done++;
    //  printf("thread %d done!\n",AC->thread_number);
    pthread_mutex_unlock(&lock);
}
void MyNLP::task_Density(void *args){
    args_class* AC = (args_class*)args;
    assert(AC != nullptr);

    //  printf("this is thread %d",AC->thread_number);
    // update with x
//    AC->args_mynlp->BoundX(AC->size_x, AC->args_mynlp->x, \
                           AC->args_mynlp->x_l, AC->args_mynlp->x_u);
//    AC->args_mynlp->UpdateExpValueForEachCell(AC->size_x, AC->args_mynlp->x,\
                                              AC->args_mynlp->_expX, AC->args_mynlp->_alpha);
//    AC->args_mynlp->UpdateExpValueForEachPin(AC->size_x, AC->args_mynlp->x,\
                                             AC->args_mynlp->_expX, AC->args_mynlp->_alpha);
//    AC->args_mynlp->UpdateNetsSumExp(AC->args_mynlp->x,\
                                     AC->args_mynlp->_expX);
//    AC->args_mynlp->UpdateBlockPosition(AC->args_mynlp->x);
//    AC->args_mynlp->UpdateDensityGrid(AC->size_x,AC->args_mynlp->x);
    AC->args_mynlp->UpdatePotentialGrid(AC->args_mynlp->x);
    SGD_Density(AC->args_mynlp,AC->size_x,AC->args_snc_Id,AC->bin,AC->stepsize,AC->maxloop);// need x back

    //vector<double> temp;
    //for(size_t i = 0;i < 2*AC->args_mynlp->m_pDB->m_modules.size();i++)
        //temp.push_back(AC->args_mynlp->x[i]);
    //save_x_threads.push_back(temp);
    pthread_mutex_lock(&lock);
    for(size_t i = 0;i < 2*AC->args_mynlp->m_pDB->m_modules.size();i++){
        Xmerge[i] += AC->args_mynlp->x[i];
        //printf("%4.3f\n",AC->args_mynlp->x[i]);
    }
    //save_x_threads.push_back(temp);
    job_done++;
    if(job_done == param.Den_thread_size + param.WL_thread_size)
        pthread_cond_signal(&cond_lock);
    if(param.debug){
        printf("weightwire = %.4f,density = %.4f %4d job(s) done,thread %4d done!\n",\
               AC->args_mynlp->_weightWire,\
               AC->args_mynlp->_weightDensity,\
               job_done,AC->thread_number);
    }
    pthread_mutex_unlock(&lock);
}
void MyNLP::task_WireLength(void *args){
    args_class* AC = (args_class*)args;
    assert(AC != nullptr);

    //  printf("this is thread %d",AC->thread_number);
    // update with x
//    AC->args_mynlp->BoundX(AC->size_x, AC->args_mynlp->x,\
                           AC->args_mynlp->x_l, AC->args_mynlp->x_u);
    AC->args_mynlp->UpdateExpValueForEachCell(AC->size_x, AC->args_mynlp->x,
                                              AC->args_mynlp->_expX, AC->args_mynlp->_alpha);
    AC->args_mynlp->UpdateExpValueForEachPin(AC->size_x, AC->args_mynlp->x,\
                                             AC->args_mynlp->_expX, AC->args_mynlp->_alpha); // not used for now
    AC->args_mynlp->UpdateNetsSumExp(AC->args_mynlp->x,
                                     AC->args_mynlp->_expX);
 //   AC->args_mynlp->UpdateBlockPosition(AC->args_mynlp->x);
 //   AC->args_mynlp->UpdateDensityGrid(AC->size_x,AC->args_mynlp->x);
//    AC->args_mynlp->UpdatePotentialGrid(AC->args_mynlp->x);


    SGD_WireLength(AC->args_mynlp,AC->size_x,AC->args_snc_Id,AC->stepsize,AC->maxloop);// need x back

    //vector<double> temp;
   // for(size_t i = 0;i < 2*AC->args_mynlp->m_pDB->m_modules.size();i++)
       // temp.push_back(AC->args_mynlp->x[i]);
    //save_x_threads.push_back(temp);
    pthread_mutex_lock(&lock);

    for(size_t i = 0;i < 2*AC->args_mynlp->m_pDB->m_modules.size();i++){
        Xmerge[i] += AC->args_mynlp->x[i];
        //printf("%4.3f\n",AC->args_mynlp->x[i]);
    }
    //save_x_threads.push_back(temp);
    job_done++;
    if(job_done == param.Den_thread_size + param.WL_thread_size)
        pthread_cond_signal(&cond_lock);
    if(param.debug){
        printf("weightwire = %.4f,density = %.4f %4d job(s) done,thread %4d done!\n",\
               AC->args_mynlp->_weightWire,\
               AC->args_mynlp->_weightDensity,\
               job_done,AC->thread_number);
    }
    pthread_mutex_unlock(&lock);
}
void shell_task_CG(void (*fp)(void*),void *args){
    /*
    *fp(args);
    vector<double> temp;
    args_class* AC = (args_class*)args;
    for(size_t i = 0;i < 2*AC->args_mynlp->m_pDB->m_nModules;i++)
        temp.push_back(AC->args_mynlp->x[i]);
    pthread_mutex_lock(&lock);
    done++;
    save_x_threads.push_back(temp);
    pthread_mutex_unlock(&lock);
    */

}
bool   MyNLP::nets_cluster_analysics(){
    double time_NCA = seconds();
    unsigned int total_nets = m_pDB->m_nets.size();
    nets_cluster.clear();
    vector<bool> is_nets_record;

    is_nets_record.resize(total_nets,false);
    for(unsigned int netId = 0;netId < total_nets;netId++){

        vector<int> nearby_nets;

        if(is_nets_record[netId]) // has already been recorded
            continue;
        int module_size = m_pDB->m_nets[netId].size();
        if(module_size == 0){// empty nets(without modules) -> pin?
            is_nets_record[netId] = true;
            continue;
        };
        for(int i = 0;i < module_size;i++){
            int pinId = m_pDB->m_nets[netId][i];
            int moduleId = m_pDB->m_pins[pinId].moduleId;
            for(unsigned int i_m_netsId = 0;i_m_netsId < m_pDB->m_modules[moduleId].m_netsId.size();i_m_netsId++){
                int m_netsId = m_pDB->m_modules[moduleId].m_netsId[i_m_netsId];
                if( (!is_nets_record[ m_netsId ]) ){
                    is_nets_record[ m_netsId ] = true;
                    nearby_nets.push_back(m_netsId);
                }
            }
        }

        nets_cluster.push_back(nearby_nets);
    }
    if(param.bShow)
        printf("nets cluster time: %f\n",seconds()-time_NCA);
    //test
    if(0) {   int counts = 0;
        for(int i = 0; i < nets_cluster.size(); i++){
            printf("cluster %4d,size = %4d\n",i,nets_cluster[i].size());
            counts += nets_cluster[i].size();
        }
        assert(counts == total_nets);
    }

    return true;
}
bool MyNLP::NCA(int moduleId, vector<bool>& is_nets_record, vector<int>& module_net){

    unsigned int nets_size = m_pDB->m_modules[moduleId].m_netsId.size();
    for(unsigned int i_nets = 0;i_nets < nets_size;i_nets++){
        int m_netsId = m_pDB->m_modules[moduleId].m_netsId[i_nets];
        if(!is_nets_record[m_netsId]){
            module_net.push_back(m_netsId);
            is_nets_record[m_netsId] = true;
        }
    }

    return true;
}
bool MyNLP::NCA_deeper(){
    double time_NCA = seconds();

    unsigned int total_nets = m_pDB->m_nets.size();
    nets_cluster.clear();
    vector<bool> is_nets_record;

    is_nets_record.resize(total_nets,false);
    for(unsigned int netId = 0;netId < total_nets;netId++){

        vector<int> nearby_nets;

        if(is_nets_record[netId]) // has already been recorded
            continue;
        int module_size = m_pDB->m_nets[netId].size();
        if(module_size == 0){// empty nets(without modules) -> pin?
            is_nets_record[netId] = true;
            //      nearby_nets.push_back(netId);
            //           nets_cluster.push_back(nearby_nets);
            continue;
        };
        for(int i = 0;i < module_size;i++){
            int pinId = m_pDB->m_nets[netId][i];
            int moduleId = m_pDB->m_pins[pinId].moduleId;
            NCA(moduleId,is_nets_record,nearby_nets);
            for(unsigned i_nets = 0;i_nets < m_pDB->m_modules[moduleId].m_netsId.size();i_nets++){
                int m_netsId = m_pDB->m_modules[moduleId].m_netsId[i_nets];
                for(unsigned i_moduleId = 0;  i_moduleId < m_pDB->m_nets[m_netsId].size();i_moduleId++){
                    int pinId2 = m_pDB->m_nets[m_netsId][i_moduleId];
                    int moduleId2 = m_pDB->m_pins[pinId2].moduleId;
                    NCA(moduleId2,is_nets_record,nearby_nets);
                }

            }
        }
        nets_cluster.push_back(nearby_nets);

        //test
    }
    if(0){
        unsigned int netcount = 0,modulecount = 0;
        for(int i = 0;i < nets_cluster.size();i++){
            netcount += nets_cluster[i].size();
            for(int j = 0;j < nets_cluster[i].size();j++){
                int netId = nets_cluster[i][j];
                modulecount += m_pDB->m_nets[netId].size();
            }
        }
        printf("net count = %d,module count = %d\n",netcount,modulecount);

        if(param.bShow)
            printf("\nnets cluster time: %f,cluster size = %d\n",
                   seconds()-time_NCA,
                   nets_cluster.size());
    }
    return true;
}

bool   MyNLP::smooth_part_eval_grad_f(const double* x,double* part_grad_f,int i_nets,int num_nets)
{
    double time_used = seconds();

    for ( unsigned int i = 0; i < m_pDB->m_modules.size(); i++ )
    {
        part_grad_wire[2 * i] = 0;
        part_grad_wire[2 * i + 1] = 0;
    }

    // grad WL
    if( _weightWire > 0 )
        for(unsigned int i = 0;i < num_nets;i++){
            int netId = i + i_nets;
            for(unsigned int j = 0; j < m_pDB->m_nets[netId].size();j++){
                int pinId = m_pDB->m_nets[netId][j];
                int moduleId = m_pDB->m_pins[pinId].moduleId;
                if( m_pDB->m_modules[moduleId].m_isFixed || m_pDB->m_modules[moduleId].m_netsId.size() == 0 )
                    continue;
                int selfPinId = netId;
                assert(m_usePin[moduleId] == false);
                if(m_usePin[moduleId]){
                    assert( selfPinId != -1 );
                    double xx = x[ 2*moduleId]   + m_pDB->m_pins[ selfPinId ].xOff;
                    double yy = x[ 2*moduleId+1 ] + m_pDB->m_pins[ selfPinId ].yOff;
                    grad_wire[ 2*moduleId ] +=
                            m_nets_sum_p_x_pos[netId]     * _expPins[2*selfPinId] / xx -
                            m_nets_sum_p_inv_x_neg[netId] / _expPins[2*selfPinId] / xx;
                    grad_wire[ 2*moduleId+1 ] +=
                            m_nets_sum_p_y_pos[netId]     * _expPins[2*selfPinId+1] / yy -
                            m_nets_sum_p_inv_y_neg[netId] / _expPins[2*selfPinId+1] / yy;

                }
                else{
                    double xx = x[ 2*moduleId ];
                    double yy = x[ 2*moduleId+1 ];
                    xx *= m_posScale;
                    yy *= m_posScale;

                    grad_wire[ 2*moduleId ] +=
                            m_nets_sum_p_x_pos[netId]     * _expX[2*moduleId] / xx  -
                            m_nets_sum_p_inv_x_neg[netId] / _expX[2*moduleId] / xx;
                    grad_wire[ 2*moduleId+1 ] +=
                            m_nets_sum_p_y_pos[netId]     * _expX[2*moduleId+1] / yy -
                            m_nets_sum_p_inv_y_neg[netId] / _expX[2*moduleId+1] / yy;
                }
            }
        }
    time_grad_wl += seconds() - time_used;

    // grad Density
    double time_start_2 = seconds();


    double partGradDensityX;
    double partGradDensityY;
    for(int j = 0;j < num_nets;j++){
        for(int i = 0;i < m_pDB->m_nets[i_nets+j].size();i++){
            int pinId    =  m_pDB->m_nets[i_nets+j][i];
            int moduleId =  m_pDB->m_pins[pinId].moduleId;

            if( m_pDB->m_modules[moduleId].m_isFixed == true)
                continue;
            GetPotentialGrad( x, moduleId, partGradDensityX, partGradDensityY );
            part_grad_potential[2 * moduleId]     = partGradDensityX;
            part_grad_potential[2 * moduleId + 1] = partGradDensityY;

        }
    }



    time_grad_potential += seconds() - time_start_2;

    // compute total fouce
    for(int j = 0,k = 0;j < num_nets;j++){
        for( int i = 0; i < m_pDB->m_nets[i_nets + j].size() ; i++ ){
            int pinId    =  m_pDB->m_nets[i_nets + j][i];
            int moduleId =  m_pDB->m_pins[pinId].moduleId;
            part_grad_f[2*k]   = _weightDensity * part_grad_potential[2*moduleId]   + part_grad_wire[2*moduleId]   * _weightWire;
            part_grad_f[2*k+1] = _weightDensity * part_grad_potential[2*moduleId+1] + part_grad_wire[2*moduleId+1] * _weightWire;
            k++;
        }
    }

    time_grad_f += seconds()-time_used;
    return true;
}
bool   MyNLP::smooth_part_eval_grad_f(const double* x,
                                      const vector<int> cur_nets_cluster_modules_bounds, const vector<int> cur_black_module,
                                      vector<int>& RegionModule,
                                      vector<double> &grad_f, int flag)
{
    for(int i = 0; i < grad_f.size();i++){
        grad_wire[i] = 0;
        grad_potential[i] = 0;
    }

    double time_used = seconds();
    // grad WL
    if( _weightWire > 0 )	//TEST
        for( unsigned int i = 0; i< cur_black_module.size(); i++ )	// for each block
        {
            int moduleId = cur_black_module[i];
            if( m_pDB->m_modules[moduleId].m_isFixed || m_pDB->m_modules[moduleId].m_netsId.size() == 0 )
                continue;

            for( unsigned int j=0; j<m_pDB->m_modules[moduleId].m_netsId.size(); j++ )
            {
                // for each net connecting to the block
                int netId = m_pDB->m_modules[moduleId].m_netsId[j];
                if( m_pDB->m_nets[netId].size() == 0 ) // floating-module
                    continue;

                // TODO: modification for LEF/DEF input
                // no floating pin for bookshelf format
                //if( m_pDB->m_nets[netId].size() == 0 )
                //	continue;

                int selfPinId = m_moduleNetPinId[moduleId][j];

                if( m_usePin[moduleId] )
                {
                    assert( selfPinId != -1 );
                    double xx = x[ 2*moduleId ]   + m_pDB->m_pins[ selfPinId ].xOff;
                    double yy = x[ 2*moduleId+1 ] + m_pDB->m_pins[ selfPinId ].yOff;
                    xx *= m_posScale;
                    yy *= m_posScale;

                    grad_wire[ 2*moduleId ] +=
                            m_nets_sum_p_x_pos[netId]     * _expPins[2*selfPinId] / xx -
                            m_nets_sum_p_inv_x_neg[netId] / _expPins[2*selfPinId] / xx;
                    grad_wire[ 2*moduleId+1 ] +=
                            m_nets_sum_p_y_pos[netId]     * _expPins[2*selfPinId+1] / yy -
                            m_nets_sum_p_inv_y_neg[netId] / _expPins[2*selfPinId+1] / yy;
                }
                else
                {
                    double xx = x[ 2*moduleId ];
                    double yy = x[ 2*moduleId+1 ];
                    xx *= m_posScale;
                    yy *= m_posScale;

                    grad_wire[ 2*moduleId ]   +=
                            m_nets_sum_p_x_pos[netId]     * _expX[2*moduleId] / xx  -
                            m_nets_sum_p_inv_x_neg[netId] / _expX[2*moduleId] / xx;
                    grad_wire[ 2*moduleId+1 ] +=
                            m_nets_sum_p_y_pos[netId]     * _expX[2*moduleId+1] / yy -
                            m_nets_sum_p_inv_y_neg[netId] / _expX[2*moduleId+1] / yy;
                }

            } // for each pin in the module

        } // for each module
    time_grad_wl += seconds() - time_used;

    // grad Density
    double time_start_2 = seconds();



    int gx1 = cur_nets_cluster_modules_bounds[0];
    int gx2 = cur_nets_cluster_modules_bounds[1];
    int gy1 = cur_nets_cluster_modules_bounds[2];
    int gy2 = cur_nets_cluster_modules_bounds[3];

    if(RegionModule.size() == 0){
        RegionModule.clear();
        int size_module = m_pDB->m_modules.size();
        for(int i = 0;i < size_module;i++){
            double partGradDensityX = 0;
            double partGradDensityY = 0;
            UpdateOneModulePotentialGrid(i,
                                         RegionModule,
                                         gx1,gx2,gy1,gy2,
                                         partGradDensityX,partGradDensityY);
            grad_potential[2*i]   = partGradDensityX;
            grad_potential[2*i+1] = partGradDensityY;
        }
    }
    else{
        int size_Rmodule = RegionModule.size();
        for(int i = 0;i < size_Rmodule;i++){
            double partGradDensityX = 0;
            double partGradDensityY = 0;
            int moduleId = RegionModule[i];
            //           UpdateOneModulePotentialGrid(moduleId,partGradDensityX,partGradDensityY);

            UpdateOneModulePotentialGrid_adjust(moduleId,
                                                gx1,gx2,gy1,gy2,
                                                partGradDensityX,partGradDensityY);
            grad_potential[2*moduleId]   = partGradDensityX;
            grad_potential[2*moduleId+1] = partGradDensityY;

        }
    }



    time_grad_potential += seconds() - time_start_2;

    // compute total fouce
    for(int i = 0;i < RegionModule.size();i++){
        int moduleId = RegionModule[i];
        grad_f[2*moduleId]   = _weightDensity*grad_potential[2*moduleId]
                + _weightWire*grad_wire[2*moduleId];
        grad_f[2*moduleId+1] = _weightDensity*grad_potential[2*moduleId+1]
                + _weightWire*grad_wire[2*moduleId+1];
    }

    time_grad_f += seconds()-time_used;
    return true;
}

void   MyNLP::part_eval_grad_f(const double* x, double* part_grad_f, int i_net, int num_nets){
    double time_used = seconds();

    for ( unsigned int i = 0; i < m_pDB->m_modules.size(); i++ )
    {
        part_grad_wire[2 * i] = 0;
        part_grad_wire[2 * i + 1] = 0;
    }

    if( _weightWire > 0 )	//TEST
        for(int i = 0;i < num_nets;i++){
            if (m_pDB->m_nets[i_net+i].size() == 0) continue;
            if (m_pDB->m_nets[i_net+i].size() == 1)
            {
                //  sg_updateSubGradWire_netBig_s(i, x);
                continue;
            }

            if (m_pDB->m_nets[i_net+i].size() == 2)
            {
                sg_updateSubGradWire_net2_s(i_net + i, x,true);
                // cnt2++;
                continue;
            }

            if (m_pDB->m_nets[i_net+i].size() == 3)
            {
                sg_updateSubGradWire_netBig_s(i_net + i, x,true);
                // cnt3++;
                continue;
            }

            if (m_pDB->m_nets[i_net+i].size() > 3)
            {
                sg_updateSubGradWire_netBig_s(i_net + i, x,true);
                //cnt4++;
                continue;
            }
        }

    time_grad_wl += seconds() - time_used;

    // grad Density
    double time_start_2 = seconds();

    double partGradDensityX;
    double partGradDensityY;
    for(int j = 0;j < num_nets;j++){
        for(int i = 0;i < m_pDB->m_nets[i_net+j].size();i++){
            int pinId    =  m_pDB->m_nets[i_net+j][i];
            int moduleId =  m_pDB->m_pins[pinId].moduleId;

            if( m_pDB->m_modules[moduleId].m_isFixed == true)
                continue;
            GetPotentialGrad( x, moduleId, partGradDensityX, partGradDensityY );
            part_grad_potential[2 * moduleId]     = partGradDensityX;
            part_grad_potential[2 * moduleId + 1] = partGradDensityY;

        }
    }


    time_grad_potential += seconds() - time_start_2;

    // compute total fouce
    for(int j = 0,k = 0;j < num_nets;j++){
        for( int i = 0; i < m_pDB->m_nets[i_net + j].size() ; i++ ){
            int pinId    =  m_pDB->m_nets[i_net + j][i];
            int moduleId =  m_pDB->m_pins[pinId].moduleId;
            part_grad_f[2*k]   = _weightDensity * part_grad_potential[2*moduleId]   + part_grad_wire[2*moduleId]   * _weightWire;
            part_grad_f[2*k+1] = _weightDensity * part_grad_potential[2*moduleId+1] + part_grad_wire[2*moduleId+1] * _weightWire;
            k++;
        }
    }

    time_grad_f += seconds()-time_used;
}
double MyNLP::part_CalcHPWL(const double* x,int i_nets)
{
    double HPWL = 0;
    double maxX = -10000;
    double maxY = -10000;
    double minX = DBL_MAX;
    double minY = DBL_MAX;
    double cx, cy;
    int pID, mID;

    for (unsigned int j = 0; j < m_pDB->m_nets[i_nets].size(); j++)
    {
        pID = m_pDB->m_nets[i_nets][j];
        mID = m_pDB->m_pins[pID].moduleId;
        assert( pID < m_pDB->m_pins.size() );
        assert( mID < m_pDB->m_modules.size() );
        cx = x[ 2 * mID ];
        cy = x[ 2 * mID + 1 ];

        minX = min( minX, cx);
        maxX = max( maxX, cx);
        minY = min( minY, cy);
        maxY = max( maxY, cy);
    }

    assert(maxX >= minX);
    assert(maxY >= minY);

    HPWL += ( (maxX - minX) + (maxY - minY) );


    return HPWL;

}
bool   MyNLP::test_part_eval_grad_f(const double* x){
    int n = 2*m_pDB->m_modules.size();
    double* sub_grad_f = new double [n];
    memset( sub_grad_f, 0, sizeof(double)*n );
    sg_eval_grad_f(n,x,sub_grad_f);

    double* part_grad_f = new double [n];
    memset( part_grad_f, 0, sizeof(double)*n );

    for(unsigned int i_nets = 0;i_nets < m_pDB->m_nets.size();i_nets++){
        int current_module_size = m_pDB->m_nets[i_nets].size();
        if(current_module_size == 0)
            continue;
        double* cur_grad_f = new double [2*current_module_size];
        memset( cur_grad_f, 0, sizeof(double)*2*current_module_size);
        part_eval_grad_f(x,part_grad_f,i_nets);

        for(int i = 0;i < current_module_size;i++){
            int pinId =  m_pDB->m_nets[i_nets][i];
            int moduleId =  m_pDB->m_pins[pinId].moduleId;
            part_grad_f[2*moduleId] += part_grad_f[2*i];
            part_grad_f[2*moduleId+1] += part_grad_f[2*i+1];
        }
        delete[] cur_grad_f;
    }
    double sum_diff = 0;
    int diff_count = 0;
    vector<int> diff_moduleID;
    printf("print sub_grad_f,part_grad_f\n");
    for(int i = 0;i < n;i++){
        printf("%d,sub_grad_f = %f ,part_grad_f = %f \n",i,sub_grad_f[i],part_grad_f[i]);
        double diff = abs( sub_grad_f[i] - part_grad_f[i]);
        sum_diff += diff;
        if( diff < 0.0001 )
            continue;
        else{
            //            printf("error: part_eval_grad_f fail!\n,diff = %f,sum_diff = %f",diff,sum_diff);
            //            return false;
            diff_moduleID.push_back(i);
            diff_count++;

        }
    }

    for(int i = 0; i < diff_count;i++){
        printf("%d,sub_grad[%d] = %f,part_grad[%d] = %f\n",i,diff_moduleID[i],
               sub_grad_f[ diff_moduleID[i] ],diff_moduleID[i],part_grad_f[ diff_moduleID[i] ]);
    }
    printf("%d\n",diff_count);
    return true;
}
void   MyNLP::OneModulePotentialGrid(const int cur_moduleId, const double* position, vector< vector<double> >& potential_record){
    if( m_pDB->m_modules[cur_moduleId].m_isOutCore )
        return;

    // preplaced blocks are stored in m_basePotential
    if( m_pDB->m_modules[cur_moduleId].m_isFixed )
        return;
    //for movable modules
    int gx, gy;
    double cellX = position[cur_moduleId*2];//  module i central position
    double cellY = position[cur_moduleId*2+1];
    double potentialRX = _potentialRX;
    double potentialRY = _potentialRY;
    double width  = m_pDB->m_modules[cur_moduleId].m_width;
    double height = m_pDB->m_modules[cur_moduleId].m_height;
    double len1 = width *0.5 + potentialRX;
    double len2 = height*0.5 + potentialRY;
    double left   = cellX - len1;
    double bottom = cellY - len2;
    double right  = cellX + len1;//right = cellX + width*0.5 + potentialRX
    double top    = cellY + len2;//top = cellY + height*0.5 + potentialRY
    if( left   < m_pDB->m_coreRgn.left )     left   = m_pDB->m_coreRgn.left;
    if( bottom < m_pDB->m_coreRgn.bottom )   bottom = m_pDB->m_coreRgn.bottom;
    if( top    > m_pDB->m_coreRgn.top )      top    = m_pDB->m_coreRgn.top;
    if( right  > m_pDB->m_coreRgn.right )    right  = m_pDB->m_coreRgn.right;
    GetClosestGrid( left, bottom, gx, gy );

    double totalPotential = 0;
    vector< potentialStruct > potentialList;
    int gxx, gyy;
    double xx, yy;

    //// TEST (convert to std-cell)
    if( height < m_potentialGridHeight && width < m_potentialGridWidth )
        width = height = 0;

    for( gxx = gx, xx = GetXGrid(gxx); xx<=right ; gxx++, xx+=m_potentialGridWidth )
    {     //xx =  m_pDB->m_coreRgn.left + gx * m_potentialGridWidth + 0.5 * m_potentialGridWidth
        for( gyy = gy, yy = GetYGrid(gyy); yy<=top ; gyy++, yy+=m_potentialGridHeight )
        { // yy = m_pDB->m_coreRgn.bottom + gy * m_potentialGridHeight + 0.5 * m_potentialGridHeight
            double potential = GetPotential( cellX, xx, potentialRX, width ) *
                    GetPotential( cellY, yy, potentialRY, height );
            if( potential > 0 )
            {
                totalPotential += potential;
                potentialList.push_back( potentialStruct( gxx, gyy, potential ) );
            }
        }
    }

    // normalize the potential so that total potential equals the cell area
    double scale = m_pDB->m_modules[cur_moduleId].m_area / totalPotential;
    //printf( "totalPotential = %f\n", totalPotential );

    _cellPotentialNorm[cur_moduleId] = scale;	    // normalization factor for the cell i
    vector< potentialStruct >::const_iterator ite;
    for( ite=potentialList.begin(); ite!=potentialList.end(); ++ite )
    {

        potential_record[ ite->gx ][ ite->gy ] += ite->potential * scale;
    }
}
void   MyNLP::part_UpdatePotentialGrid(const double *x,const double *x_bak,const vector<int>& moduleId,const int num_module)
{
    double time_start = seconds();
    int cur_moduleId = 0;

    vector< vector<double> > last_potential;
    last_potential.resize( m_potentialGridSize );
    for(int i = 0;i < m_potentialGridSize;i++)
        last_potential[i].resize(m_potentialGridSize,0);

    for( int i=0; i<(int)num_module; i++ )
    {
        // for each cell. cell ci coordinate is ( x[i*2], x[i*2+1] )
        cur_moduleId = moduleId[i];
        OneModulePotentialGrid(cur_moduleId,x_bak,last_potential  );
        OneModulePotentialGrid(cur_moduleId,x,    m_gridPotential );
    } // for each cell

    for(int i = 0;i < m_potentialGridSize;i++)
        for(int j = 0;j < m_potentialGridSize;j++){
            m_gridPotential[i][j] -= last_potential[i][j];
            //           assert( m_gridPotential[i][j] == last_potential[i][j]);
        }
    time_up_potential += seconds() - time_start;

}
void   MyNLP::part_UpdatePotentialGrid(const double *x){
    double time_start = seconds();
    ClearPotentialGrid();// set m_gridPotential to 0
    for(int i = 0;i < m_pDB->m_modules.size();i++)
        OneModulePotentialGrid(i,x,m_gridPotential);
    // for each cell
    time_up_potential += seconds() - time_start;
}

void MyNLP::getBlackModule(int nets_clusterId){

    vector<int> black_module;
    vector<bool> isrecord;
    isrecord.resize(m_pDB->m_modules.size(),false);
    for(unsigned int i = 0;i < nets_cluster[nets_clusterId].size();i++){
        int netId = nets_cluster[nets_clusterId][i];
        for(unsigned int j = 0;j < m_pDB->m_nets[netId].size();j++){
            int pinId = m_pDB->m_nets[netId][j];
            int moduleId = m_pDB->m_pins[pinId].moduleId;
            if( !isrecord[moduleId] ){
                black_module.push_back(moduleId);
                isrecord[moduleId] = true;
            }

        }
    }
    nets_black_modules.push_back(black_module);
}
void MyNLP::updateBlackBounds(int net_cluserId){
    int maxX = -10000;
    int maxY = -10000;
    int minX = 10000;
    int minY = 10000;
    nets_cluster_modules_bounds[net_cluserId].clear();
    double potentialRX = _potentialRX;// 2 * m_PotentialGridWidth;
    double potentialRY = _potentialRY;// 2 * m_PotentialGridHeight
    for(unsigned int i = 0;i < nets_cluster[net_cluserId].size();i++){
        int ncId = nets_cluster[net_cluserId][i];
        for(unsigned int k = 0;k < m_pDB->m_nets[ncId].size();k++){
            int pinId = m_pDB->m_nets[ncId][k];
            int moduleId = m_pDB->m_pins[pinId].moduleId;
            if(m_pDB->m_modules[moduleId].m_isFixed == true)
                continue;
            double cellX = x[moduleId*2];
            double cellY = x[moduleId*2+1];

#if 0
            if (0 && m_usePin[pinId] == true)  // don't know what it means
            {
                cellX = x[2 * moduleId]     + m_pDB->m_pins[pinId].xOff;
                cellY = x[2 * moduleId + 1] + m_pDB->m_pins[pinId].yOff;
            }
#endif

            double width  = m_pDB->m_modules[moduleId].m_width;
            double height = m_pDB->m_modules[moduleId].m_height;
            double left   = cellX - width * 0.5  - potentialRX;
            double bottom = cellY - height * 0.5 - potentialRY;
            double right  = cellX + (cellX - left);
            double top    = cellY + (cellY - bottom);
            if( left   < m_pDB->m_coreRgn.left )     left   = m_pDB->m_coreRgn.left;
            if( bottom < m_pDB->m_coreRgn.bottom )   bottom = m_pDB->m_coreRgn.bottom;
            if( top    > m_pDB->m_coreRgn.top )      top    = m_pDB->m_coreRgn.top;
            if( right  > m_pDB->m_coreRgn.right )    right  = m_pDB->m_coreRgn.right;
            int gx1,gx2,gy1,gy2;
            GetClosestGrid( left, bottom, gx1, gy1 );
            GetClosestGrid( right, top, gx2, gy2 );


            maxX = (maxX > gx2) ? maxX:gx2;
            maxY = (maxY > gy2) ? maxY:gy2;
            minX = (minX < gx1) ? minX:gx1;
            minY = (minY < gy1) ? minY:gy1;

        }


    }
#if 0
    assert(minX >= 0 && minX <= maxX);
    assert(minY >= 0 && minY <= maxY);
    assert(maxX <= m_gridPotential.size());
    assert(maxY <= m_gridPotential.size());
#endif
    if(maxX == m_gridPotential.size())
        maxX -=1;
    if(maxY == m_gridPotential.size())
        maxY -=1;
    nets_cluster_modules_bounds[net_cluserId].push_back(minX);
    nets_cluster_modules_bounds[net_cluserId].push_back(maxX);
    nets_cluster_modules_bounds[net_cluserId].push_back(minY);
    nets_cluster_modules_bounds[net_cluserId].push_back(maxY);
}
void MyNLP::updateBlackBounds(int net_cluserId,
                              const vector< vector<int> >& cur_nets_cluster,
                              const vector< vector<int> >& cur_nets_black_modules,
                              vector< vector<int> >& cur_nets_cluster_modules_bounds){
    int maxX = -10000;
    int maxY = -10000;
    int minX = 10000;
    int minY = 10000;
    cur_nets_cluster_modules_bounds[net_cluserId].clear();
    double potentialRX = _potentialRX;// 2 * m_PotentialGridWidth;
    double potentialRY = _potentialRY;// 2 * m_PotentialGridHeight
    for(unsigned int i = 0;i < cur_nets_cluster[net_cluserId].size();i++){
        int ncId = cur_nets_cluster[net_cluserId][i];
        for(unsigned int k = 0;k < m_pDB->m_nets[ncId].size();k++){
            int pinId = m_pDB->m_nets[ncId][k];
            int moduleId = m_pDB->m_pins[pinId].moduleId;
            if(m_pDB->m_modules[moduleId].m_isFixed == true)
                continue;
            double cellX = x[moduleId*2];
            double cellY = x[moduleId*2+1];

#if 0
            if (0 && m_usePin[pinId] == true)  // don't know what it means
            {
                cellX = x[2 * moduleId]     + m_pDB->m_pins[pinId].xOff;
                cellY = x[2 * moduleId + 1] + m_pDB->m_pins[pinId].yOff;
            }
#endif

            double width  = m_pDB->m_modules[moduleId].m_width;
            double height = m_pDB->m_modules[moduleId].m_height;
            double left   = cellX - width * 0.5  - potentialRX;
            double bottom = cellY - height * 0.5 - potentialRY;
            double right  = cellX + (cellX - left);
            double top    = cellY + (cellY - bottom);
            if( left   < m_pDB->m_coreRgn.left )     left   = m_pDB->m_coreRgn.left;
            if( bottom < m_pDB->m_coreRgn.bottom )   bottom = m_pDB->m_coreRgn.bottom;
            if( top    > m_pDB->m_coreRgn.top )      top    = m_pDB->m_coreRgn.top;
            if( right  > m_pDB->m_coreRgn.right )    right  = m_pDB->m_coreRgn.right;
            int gx1,gx2,gy1,gy2;
            GetClosestGrid( left, bottom, gx1, gy1 );
            GetClosestGrid( right, top, gx2, gy2 );


            maxX = (maxX > gx2) ? maxX:gx2;
            maxY = (maxY > gy2) ? maxY:gy2;
            minX = (minX < gx1) ? minX:gx1;
            minY = (minY < gy1) ? minY:gy1;

        }


    }
#if 0
    assert(minX >= 0 && minX <= maxX);
    assert(minY >= 0 && minY <= maxY);
    assert(maxX <= m_gridPotential.size());
    assert(maxY <= m_gridPotential.size());
#endif
    if(maxX == m_gridPotential.size())
        maxX -=1;
    if(maxY == m_gridPotential.size())
        maxY -=1;
    cur_nets_cluster_modules_bounds[net_cluserId].push_back(minX);
    cur_nets_cluster_modules_bounds[net_cluserId].push_back(maxX);
    cur_nets_cluster_modules_bounds[net_cluserId].push_back(minY);
    cur_nets_cluster_modules_bounds[net_cluserId].push_back(maxY);
}

void MyNLP::part_eval_grad_f(const double* x, int net_cluster_id,
                             vector<int>& RegionModule,vector<double> &grad_f,int flag){
    double part_eval_grad_f_time = seconds();
    part_grad_wire.clear();
    part_grad_potential.clear();
    grad_f.clear();
    part_grad_wire.resize(2*m_pDB->m_modules.size(),0.0);
    part_grad_potential.resize(2*m_pDB->m_modules.size(),0.0);
    grad_f.resize(2*m_pDB->m_modules.size(),0.0);


    int gx1,gx2,gy1,gy2;// get the bin we concern
    gx1 = nets_cluster_modules_bounds[net_cluster_id][0];
    gx2 = nets_cluster_modules_bounds[net_cluster_id][1];
    gy1 = nets_cluster_modules_bounds[net_cluster_id][2];
    gy2 = nets_cluster_modules_bounds[net_cluster_id][3];
    double time_used = seconds();

    if( _weightWire > 0 )	//compute the HPWL grad
        for(unsigned int i = 0;i < nets_cluster[net_cluster_id].size();i++){
            int netsId = nets_cluster[net_cluster_id][i];
            if (m_pDB->m_nets[netsId].size() == 0) continue;
            if (m_pDB->m_nets[netsId].size() == 1)
            {
                //  sg_updateSubGradWire_netBig_s(i, x);
                continue;
            }

            if (m_pDB->m_nets[netsId].size() == 2)
            {
                sg_updateSubGradWire_net2_s(netsId, x,true);
                // cnt2++;
                continue;
            }

            if (m_pDB->m_nets[netsId].size() == 3)
            {
                sg_updateSubGradWire_netBig_s(netsId, x,true);
                // cnt3++;
                continue;
            }

            if (m_pDB->m_nets[netsId].size() > 3)
            {
                sg_updateSubGradWire_netBig_s(netsId, x,true);
                //cnt4++;
                continue;
            }
        }

    time_grad_wl += seconds() - time_used;
    if(flag){
        for(unsigned int i = 0; i < m_pDB->m_modules.size();i++){//compute the density
            double partGradDensityX,partGradDensityY;
            UpdateOneModulePotentialGrid(i,RegionModule,gx1,gx2,gy1,gy2,partGradDensityX,partGradDensityY);
            part_grad_potential[2*i] = partGradDensityX;
            part_grad_potential[2*i] = partGradDensityY;
        }
    }
    else{
        for(unsigned int i = 0; i < RegionModule.size();i++){//compute the density
            double partGradDensityX,partGradDensityY;
            int moduleId = RegionModule[i];
            UpdateOneModulePotentialGrid_adjust(moduleId,
                                                // RegionModule,
                                                gx1,gx2,gy1,gy2,
                                                //  i,
                                                partGradDensityX,partGradDensityY);
            part_grad_potential[2*moduleId] = partGradDensityX;
            part_grad_potential[2*moduleId] = partGradDensityY;
        }
    }
    for(unsigned int i = 0; i < m_pDB->m_modules.size();i++){
        double temp1 = _weightWire*part_grad_wire[2*i]   + _weightDensity*part_grad_potential[2*i];
        double temp2 = _weightWire*part_grad_wire[2*i+1] + _weightDensity*part_grad_potential[2*i+1];
        grad_f[2*i]   = temp1;
        grad_f[2*i+1] = temp2;
    }

    if(0)
    {
        if(flag)
            printf("\n all:part_eval_grad_f time = %f\n",seconds() - part_eval_grad_f_time);
        else
            printf("\npart:part_eval_grad_f time = %f\n",seconds() - part_eval_grad_f_time);

    }

}
void MyNLP::UpdateOneModulePotentialGrid(int moduleId,vector<int>& RegionModule,int region_gx1,int region_gx2,int region_gy1,int region_gy2,
                                         double& partGradDensityX,double& partGradDensityY){
    // for each cell. cell ci coordinate is ( x[i*2], x[i*2+1] )
    //    param.count_updatepotential++;
    //  double time1 = seconds();
    partGradDensityX = 0;
    partGradDensityY = 0;
    if( m_pDB->m_modules[moduleId].m_isOutCore )
        return;

    // preplaced blocks are stored in m_basePotential
    if( m_pDB->m_modules[moduleId].m_isFixed )
        return;

    int gx1, gy1,gx2,gy2;
    double cellX = x[moduleId*2];
    double cellY = x[moduleId*2+1];
    double potentialRX = _potentialRX;// 2 * m_PotentialGridWidth;
    double potentialRY = _potentialRY;// 2 * m_PotentialGridHeight
    double width  = m_pDB->m_modules[moduleId].m_width;
    double height = m_pDB->m_modules[moduleId].m_height;
    double len1 = width  * 0.5 + potentialRX;
    double len2 = height * 0.5 + potentialRY;
    double left   = cellX - len1;
    double bottom = cellY - len2;
    double right  = cellX + len1;// hdgao changed
    double top    = cellY + len2;// hdgao changed
    if( left   < m_pDB->m_coreRgn.left )     left   = m_pDB->m_coreRgn.left;
    if( bottom < m_pDB->m_coreRgn.bottom )   bottom = m_pDB->m_coreRgn.bottom;
    if( top    > m_pDB->m_coreRgn.top )      top    = m_pDB->m_coreRgn.top;
    if( right  > m_pDB->m_coreRgn.right )    right  = m_pDB->m_coreRgn.right;
    GetClosestGrid( left, bottom, gx1, gy1 );
    GetClosestGrid( right,top,    gx2, gy2 );

    if( !(    ( gx1 <= region_gx2 ) && ( gy1 <= region_gy2 )
              && ( gx2 >= region_gx1 ) && ( gy2 >= region_gy1 )  ) )// hdgao changed
        return;
    RegionModule.push_back(moduleId);
    //// TEST (convert to std-cell)
    if( height < m_potentialGridHeight && width < m_potentialGridWidth )
        width = height = 0;
    int gxx_start = (gx1 < region_gx1) ? region_gx1 : gx1; //max
    int gxx_end   = (gx2 < region_gx2) ? gx2 : region_gx2; //min
    int gyy_start = (gy1 < region_gy1) ? region_gy1 : gy1; //max
    int gyy_end   = (gy2 < region_gy2) ? gy2 : region_gy2; //min

    int gxx, gyy;
    double xx, yy;

    // xx,yy the center position of bin(gx,gy)
    for(gxx = gxx_start,xx = GetXGrid(gxx); gxx <= gxx_end; gxx++, xx+=m_potentialGridWidth )
    {// m_pDB->m_coreRgn.left + gx * m_potentialGridWidth + 0.5 * m_potentialGridWidth;
        for( gyy = gyy_start,yy = GetYGrid(gyy); gyy<=gyy_end ; gyy++, yy+=m_potentialGridHeight )
        {
            double gX = 0;
            double gY = 0;
            double common = ( m_gridPotential[gxx][gyy] - m_expBinPotential[gxx][gyy] ) *_cellPotentialNorm[moduleId];
            gX = common *
                    GetGradPotential( cellX, xx, _potentialRX, width ) *
                    GetPotential(     cellY, yy, _potentialRY, height );
            gY = common *
                    GetPotential(     cellX, xx, _potentialRX, width  ) *
                    GetGradPotential( cellY, yy, _potentialRY, height );

            partGradDensityX += gX;
            partGradDensityY += gY;
        }
    }
    //    param.time_updatepotential += seconds() - time1;
}
void MyNLP::UpdateOneModulePotentialGrid_adjust( int moduleId,
                                                 //                                                 vector<int>& RegionModule,
                                                 int region_gx1,
                                                 int region_gx2,
                                                 int region_gy1,
                                                 int region_gy2,
                                                 //                                                int n_RegionModule,
                                                 double& partGradDensityX,
                                                 double& partGradDensityY){
    //    param.count_updatepotential_adjust++;
    //    double time1 = seconds();
    // for each cell. cell ci coordinate is ( x[i*2], x[i*2+1] )
    partGradDensityX = 0;
    partGradDensityY = 0;

    if( m_pDB->m_modules[moduleId].m_isOutCore )
        return ;

    // preplaced blocks are stored in m_basePotential
    if( m_pDB->m_modules[moduleId].m_isFixed )
        return ;

    int gx1, gy1,gx2,gy2;
    double cellX = x[moduleId*2];
    double cellY = x[moduleId*2+1];
    double potentialRX = _potentialRX;// 2 * m_PotentialGridWidth;
    double potentialRY = _potentialRY;// 2 * m_PotentialGridHeight
    double width  = m_pDB->m_modules[moduleId].m_width;
    double height = m_pDB->m_modules[moduleId].m_height;
    double len1 = width  * 0.5 + potentialRX;
    double len2 = height * 0.5 + potentialRY;
    double left   = cellX - len1;
    double bottom = cellY - len2;
    double right  = cellX + len1;// hdgao changed
    double top    = cellY + len2;// hdgao changed
    if( left   < m_pDB->m_coreRgn.left )     left   = m_pDB->m_coreRgn.left;
    if( bottom < m_pDB->m_coreRgn.bottom )   bottom = m_pDB->m_coreRgn.bottom;
    if( top    > m_pDB->m_coreRgn.top )      top    = m_pDB->m_coreRgn.top;
    if( right  > m_pDB->m_coreRgn.right )    right  = m_pDB->m_coreRgn.right;
    GetClosestGrid( left, bottom, gx1, gy1 );
    GetClosestGrid( right,top,    gx2, gy2 );
    if( !(    ( gx1 <= region_gx2 ) && ( gy1 <= region_gy2 )
              && ( gx2 >= region_gx1 ) && ( gy2 >= region_gy1 )  ) ){// hdgao changed{
        return ;
    }



    //// TEST (convert to std-cell)
    if( height < m_potentialGridHeight && width < m_potentialGridWidth )
        width = height = 0;
    int gxx_start = (gx1 < region_gx1) ? region_gx1 : gx1; //max
    int gxx_end   = (gx2 < region_gx2) ? gx2 : region_gx2; //min
    int gyy_start = (gy1 < region_gy1) ? region_gy1 : gy1; //max
    int gyy_end   = (gy2 < region_gy2) ? gy2 : region_gy2; //min

    int gxx,gyy;
    double xx, yy;

    for(gxx = gxx_start,xx = GetXGrid(gxx); gxx <= gxx_end; gxx++, xx+=m_potentialGridWidth )
    {// m_pDB->m_coreRgn.left + gx * m_potentialGridWidth + 0.5 * m_potentialGridWidth;
        for( gyy = gyy_start,yy = GetYGrid(gyy); gyy<=gyy_end ; gyy++, yy+=m_potentialGridHeight )
        {
            double gX = 0;
            double gY = 0;
            double common = ( m_gridPotential[gxx][gyy] - m_expBinPotential[gxx][gyy] ) *_cellPotentialNorm[moduleId];
            gX = common *
                    GetGradPotential( cellX, xx, _potentialRX, width ) *
                    GetPotential(     cellY, yy, _potentialRY, height );
            gY = common *
                    GetPotential(     cellX, xx, _potentialRX, width  ) *
                    GetGradPotential( cellY, yy, _potentialRY, height );

            partGradDensityX += gX;
            partGradDensityY += gY;
        }
    }
    //    param.time_updatepotential_adjust += seconds() - time1;
    return ;

}

void MyNLP::UpdateOneModulePotentialGrid( int moduleId,double& partGradDensityX,double& partGradDensityY){
    //    param.count_updatepotential_adjust++;
    //  double time1 = seconds();
    // for each cell. cell ci coordinate is ( x[i*2], x[i*2+1] )
    partGradDensityX = 0;
    partGradDensityY = 0;

    if( m_pDB->m_modules[moduleId].m_isOutCore )
        return ;

    // preplaced blocks are stored in m_basePotential
    if( m_pDB->m_modules[moduleId].m_isFixed )
        return ;

    int gx1, gy1,gx2,gy2;
    double cellX = x[moduleId*2];
    double cellY = x[moduleId*2+1];
    double potentialRX = _potentialRX;// 2 * m_PotentialGridWidth;
    double potentialRY = _potentialRY;// 2 * m_PotentialGridHeight
    double width  = m_pDB->m_modules[moduleId].m_width;
    double height = m_pDB->m_modules[moduleId].m_height;
    double len1 = width  * 0.5 + potentialRX;
    double len2 = height * 0.5 + potentialRY;
    double left   = cellX - len1;
    double bottom = cellY - len2;
    double right  = cellX + len1;// hdgao changed
    double top    = cellY + len2;// hdgao changed
    if( left   < m_pDB->m_coreRgn.left )     left   = m_pDB->m_coreRgn.left;
    if( bottom < m_pDB->m_coreRgn.bottom )   bottom = m_pDB->m_coreRgn.bottom;
    if( top    > m_pDB->m_coreRgn.top )      top    = m_pDB->m_coreRgn.top;
    if( right  > m_pDB->m_coreRgn.right )    right  = m_pDB->m_coreRgn.right;
    GetClosestGrid( left, bottom, gx1, gy1 );
    GetClosestGrid( right,top,    gx2, gy2 );



    //// TEST (convert to std-cell)
    if( height < m_potentialGridHeight && width < m_potentialGridWidth )
        width = height = 0;
    int gxx,gyy;
    double xx, yy;

    for(gxx = gx1,xx = GetXGrid(gxx); gxx < gx2; gxx++, xx+=m_potentialGridWidth )
    {// m_pDB->m_coreRgn.left + gx * m_potentialGridWidth + 0.5 * m_potentialGridWidth;
        for( gyy = gy1,yy = GetYGrid(gyy); gyy<  gy2 ; gyy++, yy+=m_potentialGridHeight )
        {
            double gX = 0;
            double gY = 0;
            double common = ( m_gridPotential[gxx][gyy] - m_expBinPotential[gxx][gyy] ) *_cellPotentialNorm[moduleId];
            gX = common *
                    GetGradPotential( cellX, xx, _potentialRX, width ) *
                    GetPotential(     cellY, yy, _potentialRY, height );
            gY = common *
                    GetPotential(     cellX, xx, _potentialRX, width  ) *
                    GetGradPotential( cellY, yy, _potentialRY, height );

            partGradDensityX += gX;
            partGradDensityY += gY;
        }
    }
    //   param.time_updatepotential_adjust += seconds() - time1;
    return ;

}

void MyNLP::PartAdjustForce(vector<double> &grad_f, const vector<int>& RegionModule){

    double totalGrad = 0;
    int modulesize = RegionModule.size();
    for(int i = 0 ; i < modulesize;i++){
        int moduleId = RegionModule[i];
        double value = grad_f[2*moduleId] * grad_f[2*moduleId] + grad_f[2*moduleId+1] * grad_f[2*moduleId+1];
        totalGrad += value;
    }
    double avgGrad = sqrt(totalGrad / modulesize);

    // Do truncation
    double expMaxGrad = avgGrad * truncationFactor;	// x + y
    double expMaxGradSquare = expMaxGrad * expMaxGrad;
    for( int i=0; i<modulesize; i++ )
    {
        int moduleId =  RegionModule[i];
        double valueSquare = ( grad_f[2*moduleId]*grad_f[2*moduleId] + grad_f[2*moduleId+1]*grad_f[2*moduleId+1] );
        if( valueSquare > expMaxGradSquare )
        {
            double value = sqrt( valueSquare );
            grad_f[2*moduleId]   = grad_f[2*moduleId]   * expMaxGrad / value;
            grad_f[2*moduleId+1] = grad_f[2*moduleId+1] * expMaxGrad / value;
        }
    }
}
void MyNLP::PartAdjustForce_2D(vector<double> &grad_f, const vector<int>& RegionModule){
    double totalGrad1 = 0,totalGrad2 = 0;
    int modulesize = RegionModule.size();
    for(int i = 0 ; i < modulesize;i++){
        int moduleId = RegionModule[i];
        double value1 = grad_f[2*moduleId] * grad_f[2*moduleId] ;
        double value2 = grad_f[2*moduleId+1] * grad_f[2*moduleId+1];
        totalGrad1 += value1;
        totalGrad2 += value2;
    }
    double avgGrad1 = sqrt(totalGrad1 / modulesize);
    double avgGrad2 = sqrt(totalGrad2 / modulesize);
    // Do truncation
    double expMaxGrad1 = avgGrad1 * truncationFactor;	// x + y
    double expMaxGrad2 = avgGrad2 * truncationFactor;	// x + y
    double expMaxGradSquare1 = expMaxGrad1 * expMaxGrad1;
    double expMaxGradSquare2 = expMaxGrad2 * expMaxGrad2;
    for( int i=0; i<modulesize; i++ )
    {
        int moduleId =  RegionModule[i];
        double valueSquare1 = ( grad_f[2*moduleId]  * grad_f[2*moduleId]);
        double valueSquare2 = ( grad_f[2*moduleId+1]* grad_f[2*moduleId+1] );
        if( valueSquare1 > expMaxGradSquare1 )
        {
            double value1 = sqrt( valueSquare1 );
            grad_f[2*moduleId]   = grad_f[2*moduleId]   * expMaxGrad1 / value1;
        }
        if( valueSquare2 > expMaxGradSquare2 )
        {
            double value2 = sqrt( valueSquare2 );
            grad_f[2*moduleId+1] = grad_f[2*moduleId+1] * expMaxGrad2 / value2;
        }
    }
}

void MyNLP::PartFindBeta(const vector<double>& grad_f,const vector<double>& last_grad_f,
                         const vector<int>& RegionModule,double& beta){
    // Polak-Ribiere foumula from APlace journal paper
    // NOTE:
    //   g_{k-1} = -last_grad_f
    //   g_k     = grad_f

    double l2norm = 0;
    int size = RegionModule.size();
    for( int i=0; i<size; i++ ){
        int moduleId = RegionModule[i];
        l2norm += last_grad_f[2*moduleId] * last_grad_f[2*moduleId];
        l2norm += last_grad_f[2*moduleId+1] * last_grad_f[2*moduleId+1];
    }

    double product = 0;
    for( int i=0; i<size; i++ ){
        int moduleId = RegionModule[i];
        product += grad_f[2*moduleId] * ( grad_f[2*moduleId] + last_grad_f[2*moduleId] );	// g_k^T ( g_k - g_{k-1} )
        product += grad_f[2*moduleId+1] * ( grad_f[2*moduleId+1] + last_grad_f[2*moduleId+1]);
    }
    beta = product / l2norm;
}
void MyNLP::PartFindBeta_2D(const vector<double>& grad_f,const vector<double>& last_grad_f,
                            const vector<int>& RegionModule,double& beta1,double& beta2){
    // Polak-Ribiere foumula from APlace journal paper
    // NOTE:
    //   g_{k-1} = -last_grad_f
    //   g_k     = grad_f

    double l2norm1 = 0,l2norm2 = 0;
    int size = RegionModule.size();
    for( int i=0; i<size; i++ ){
        int moduleId = RegionModule[i];
        l2norm1 += last_grad_f[2*moduleId] * last_grad_f[2*moduleId];
        l2norm2 += last_grad_f[2*moduleId+1] * last_grad_f[2*moduleId+1];
    }

    double product1 = 0,product2 = 0;
    for( int i=0; i<size; i++ ){
        int moduleId = RegionModule[i];
        product1 += grad_f[2*moduleId] * ( grad_f[2*moduleId] + last_grad_f[2*moduleId] );	// g_k^T ( g_k - g_{k-1} )
        product2 += grad_f[2*moduleId+1] * ( grad_f[2*moduleId+1] + last_grad_f[2*moduleId+1]);
    }
    beta1 = product1 / l2norm1;
    beta2 = product2 / l2norm2;
}
void MyNLP::PartLineSearch(const vector<double> &grad_f, const vector<int>& RegionModule, double& stepsize){
    int size = RegionModule.size();
    double totalGrad = 0;
    for( int i=0; i<size; i++ ){
        int moduleId = RegionModule[i];
        totalGrad += grad_f[2*moduleId] * grad_f[2*moduleId] ;
        totalGrad += grad_f[2*moduleId + 1] * grad_f[2*moduleId+1];
    }
    double avgGrad = sqrt( totalGrad / size );
    stepsize = ( m_potentialGridWidth / avgGrad ) * m_currentStep;

    return;
}
void MyNLP::PartLineSearch_2D(const vector<double> &grad_f, const vector<int>& RegionModule, double& stepsize1,double& stepsize2){
    int size = RegionModule.size();
    double totalGrad1 = 0,totalGrad2 = 0;
    for( int i=0; i<size; i++ ){
        int moduleId = RegionModule[i];
        totalGrad1 += grad_f[2*moduleId] * grad_f[2*moduleId] ;
        totalGrad2 += grad_f[2*moduleId + 1] * grad_f[2*moduleId+1];
    }
    double avgGrad1 = sqrt( totalGrad1 / size );
    double avgGrad2 = sqrt( totalGrad2 / size );
    stepsize1 = ( m_potentialGridWidth / avgGrad1 ) * m_currentStep;
    stepsize2 = ( m_potentialGridHeight / avgGrad2 ) * m_currentStep;

    return;
}
std::vector<int> MyNLP::randVector(size_t num){
    std::vector<int> result;
    result.clear();
    result.reserve(num);
    srand( (unsigned int)time(0) );
    for(int i = 0; i < num;i++)
        result.push_back(i);
    int p1,p2,temp;
    while(--num){
        p1 = num;
        p2 = rand()%num;
        temp = result[p1];
        result[p1] = result[p2];
        result[p2] = temp;
    }
    return result;
}
void MyNLP::get_center_position(double* x,double& c_x, double& c_y, const size_t netId, int flag){
    c_x = 0;
    c_y = 0;
    double total_net_area = 0;
    if(flag == 0){//ortho
        for(unsigned int i = 0; i < m_pDB->m_nets[netId].size();i++){
            int pinId    = m_pDB->m_nets[netId][i];
            int moduleId = m_pDB->m_pins[pinId].moduleId;
            c_x += x[2*moduleId]   * m_pDB->m_modules[moduleId].m_area;
            c_y += x[2*moduleId+1] * m_pDB->m_modules[moduleId].m_area;
            total_net_area += m_pDB->m_modules[moduleId].m_area;
        }
        c_x /= total_net_area;
        c_y /= total_net_area;
    }
    if(flag == 1){//mean
        for(unsigned int i = 0; i < m_pDB->m_nets[netId].size();i++){
            int pinId    = m_pDB->m_nets[netId][i];
            int moduleId = m_pDB->m_pins[pinId].moduleId;
            c_x += x[2*moduleId];
            c_y += x[2*moduleId+1];
            total_net_area += m_pDB->m_modules[moduleId].m_area;
        }
        c_x /= m_pDB->m_nets[netId].size();
        c_y /= m_pDB->m_nets[netId].size();
    }
}
void MyNLP::dummy_task(void *arg) {
    int* p = (int*)arg;
    printf("this is thread %d\n",*p);
    fflush(stdout);
    pthread_mutex_lock(&lock);
    job_done++;
    pthread_mutex_unlock(&lock);
   // usleep(10000);

}
void MyNLP::test(threadpool_t *pool){
    pool = nullptr;
}
/*===========================================make the change=================================================*/
// Density
void MyNLP::DensityGrad(vector<double>& grad,const vector<int> netsID,
                        const vector<int>& bin){
#if 0
    for(size_t i = 0;i < m_pDB->m_modules.size();i++){
        if( m_pDB->m_modules[i].m_isFixed)
            continue;
        double gradx,grady;
        GetDensityGrad(i,bin,gradx,grady);
        grad[2*i]   = gradx*_weightDensity;//*_weightDensity;
        grad[2*i+1] = grady*_weightDensity;//*_weightDensity;
    }
#endif
#if 1
    for(size_t i = 0;i < grad.size();i++)
            grad[i] = 0;
    for(size_t i = 0;i < netsID.size();i++){
        size_t netId = netsID[i];
        size_t moduleSize = m_pDB->m_nets[netId].size();
        for(size_t k = 0;k < moduleSize;k++){
            size_t pinId = m_pDB->m_nets[netId][k];
            size_t moduleId = m_pDB->m_pins[pinId].moduleId;
            double gradx,grady;
           // GetPotentialGrad(x,moduleId,gradx,grady);
            if( m_pDB->m_modules[moduleId].m_isFixed)
                continue;
            GetDensityGrad(moduleId,bin,gradx,grady);
            grad[2*moduleId]   = gradx*_weightDensity;//*_weightDensity;
            grad[2*moduleId+1] = grady*_weightDensity;//*_weightDensity;
        }
    }
#endif
}
void MyNLP::DensityGradByModule(vector<double>& grad,const vector<int>& modules,const vector<int>& bin){
    size_t modulesize = modules.size();
    for(size_t i =0;i < grad.size();i++)
        grad[i] = 0;
    for(size_t i = 0;i < modulesize;i++){
        size_t moduleId = modules[i];
        double gradx,grady;
       // GetPotentialGrad(x,moduleId,gradx,grady);
        if( m_pDB->m_modules[moduleId].m_isFixed)
            continue;
        GetDensityGrad(moduleId,bin,gradx,grady);
        grad[2*moduleId]   = gradx*_weightDensity;//*_weightDensity;
        grad[2*moduleId+1] = grady*_weightDensity;//*_weightDensity;
    }
}

void MyNLP::DensityAdjustForce(size_t size_x,vector<double>& grad,const vector<int> netsID){
    // seems conflict with SGD
    double totalGrad = 0;
    int size = size_x/2;
    for( int i=0; i<size; i++ )
    {
        double value = grad[2*i] * grad[2*i] + grad[2*i+1] * grad[2*i+1];
        totalGrad += value;
    }
    double avgGrad = sqrt( totalGrad / size );

    // Do truncation
    double expMaxGrad = avgGrad * truncationFactor;	// x + y
    double expMaxGradSquare = expMaxGrad * expMaxGrad;
  //  printf("avgGrad = %4.4f,expMaxGrad = %4.4f,expMaxGradSquare = %4.4f\n",\
           avgGrad,expMaxGrad,expMaxGradSquare);
 //   fflush(stdout);
    for( int i=0; i<size; i++ )
    {

        double valueSquare = ( grad[2*i]*grad[2*i] + grad[2*i+1]*grad[2*i+1] );
        if( valueSquare > expMaxGradSquare )
        {
            double value = sqrt( valueSquare );
            assert(value == value);
            grad[2*i]   = grad[2*i]   * expMaxGrad / value;
            grad[2*i+1] = grad[2*i+1] * expMaxGrad / value;
        }

    }
}

void MyNLP::DensityBoundX(size_t size_x){//
    for( size_t i=0; i<size_x; i++ )
    {
        if( x[i] < x_l[i] )
            x[i] = x_l[i];
        else if( x[i] > x_u[i] )
            x[i] = x_u[i];
    }
}

void MyNLP::DensityUpdatePotentialGrid(const vector<int> netsID,vector<int> bin){
    for(size_t i = 0;i < netsID.size();i++){
        size_t netId = netsID[i];
        size_t moduleSize = m_pDB->m_nets[netId].size();
        for(size_t k = 0;k < moduleSize;k++){
            size_t pinId = m_pDB->m_nets[netId][k];
            size_t moduleId = m_pDB->m_pins[pinId].moduleId;
            DensityUpdateOneModulePotentialGrid(moduleId,bin);
        }
    }
}



// wirelength
void MyNLP::WireLengthGrad(vector<double>& grad,const vector<int> netsID){

    if(_weightWire == 0)
        return;

    for(size_t i = 0;i < grad.size();i++)
            grad[i] = 0;
    for(size_t i = 0;i < netsID.size();i++){

        int netId = netsID[i];
        size_t moduleSize = m_pDB->m_nets[netId].size();
        for(size_t k = 0;k < moduleSize;k++){

            size_t pinId = m_pDB->m_nets[netId][k];
            size_t moduleId = m_pDB->m_pins[pinId].moduleId;

            if( m_pDB->m_modules[moduleId].m_isFixed || m_pDB->m_modules[moduleId].m_netsId.size() == 0 )
                continue;

            int curnetId = netId;
            if( m_pDB->m_nets[curnetId].size() == 0 ) // floating-module
                continue;

            // TODO: modification for LEF/DEF input
            // no floating pin for bookshelf format
            //if( m_pDB->m_nets[netId].size() == 0 )
            //	continue;

            int selfPinId = pinId;//m_moduleNetPinId[moduleId][j];
            assert(m_usePin[moduleId] == false);
            if( m_usePin[moduleId] ){// not used
                assert( selfPinId != -1 );
                double xx = x[ 2*moduleId ]   + m_pDB->m_pins[ selfPinId ].xOff;
                double yy = x[ 2*moduleId+1 ] + m_pDB->m_pins[ selfPinId ].yOff;
                xx *= m_posScale;
                yy *= m_posScale;

                grad[ 2*moduleId ] +=
                        m_nets_sum_p_x_pos[curnetId]     * _expPins[2*selfPinId] / xx -
                        m_nets_sum_p_inv_x_neg[curnetId] / _expPins[2*selfPinId] / xx;
                grad[ 2*moduleId+1 ] +=
                        m_nets_sum_p_y_pos[curnetId]     * _expPins[2*selfPinId+1] / yy -
                        m_nets_sum_p_inv_y_neg[curnetId] / _expPins[2*selfPinId+1] / yy;
            }
            else{
                double xx = x[ 2*moduleId ];
                double yy = x[ 2*moduleId+1 ];
                xx *= m_posScale;
                yy *= m_posScale;

                grad[ 2*moduleId ] +=
                        (m_nets_sum_p_x_pos[curnetId]     * _expX[2*moduleId] / xx  -
                        m_nets_sum_p_inv_x_neg[curnetId] / _expX[2*moduleId] / xx);
                grad[ 2*moduleId+1 ] +=
                        (m_nets_sum_p_y_pos[curnetId]     * _expX[2*moduleId+1] / yy -
                        m_nets_sum_p_inv_y_neg[curnetId] / _expX[2*moduleId+1] / yy);
            }

        }// end k

    }// end i
    for(size_t i = 0;i < grad.size();i++)
            grad[i] *= _weightWire;
}

void MyNLP::WireLengthAdjustForce(size_t size_x,vector<double>& grad,const vector<int> netsID){
    // seems conflict with SGD
    double totalGrad = 0;
    int size = size_x/2;
    for( int i=0; i<size; i++ )
    {
        double value = grad[2*i] * grad[2*i] + grad[2*i+1] * grad[2*i+1];
        totalGrad += value;
    }
    double avgGrad = sqrt( totalGrad / size );

    // Do truncation
    double expMaxGrad = avgGrad * truncationFactor;	// x + y
    double expMaxGradSquare = expMaxGrad * expMaxGrad;
    for( int i=0; i<size; i++ )
    {
        double valueSquare = ( grad[2*i]*grad[2*i] + grad[2*i+1]*grad[2*i+1] );
        if( valueSquare > expMaxGradSquare )
        {
            double value = sqrt( valueSquare );
            grad[2*i]   = grad[2*i]   * expMaxGrad / value;
            grad[2*i+1] = grad[2*i+1] * expMaxGrad / value;
        }
    }
}

void MyNLP::WireLengthBoundX(size_t size_x){
    for( size_t i=0; i<size_x; i++ ){
        if( x[i] < x_l[i] )
            x[i] = x_l[i];
        else if( x[i] > x_u[i] )
            x[i] = x_u[i];
    }
}

void MyNLP::WireLengthUpdateExpValueForCell(const vector<int> netsID, size_t size_x){
    vector<bool> isread;
    isread.resize(size_x/2,false);
    for(size_t i = 0;i < netsID.size();i++){
        size_t netId = netsID[i];
        size_t moduleSize = m_pDB->m_nets[netId].size();
        for(size_t k = 0;k < moduleSize;k++){
            size_t pinId = m_pDB->m_nets[netId][k];
            size_t moduleId = m_pDB->m_pins[pinId].moduleId;
            if(!isread[moduleId]){// don't compute twice
                isread[moduleId] = true;
                _expX[2*moduleId]   = pow( x[2*moduleId] * m_posScale, _alpha );
                _expX[2*moduleId+1] = pow( x[2*moduleId+1] * m_posScale, _alpha );
            }
        }
    }
}
void MyNLP::WireLengthUpdateExpValueForCellByModule(const vector<int> modules,size_t size_x){
    size_t moduleSize = modules.size();
    for(size_t i = 0;i < moduleSize;i++){
        size_t moduleId = modules[i];
        _expX[2*moduleId]   = pow( x[2*moduleId] * m_posScale, _alpha );
        _expX[2*moduleId+1] = pow( x[2*moduleId+1] * m_posScale, _alpha );
    }
}
void MyNLP::WireLengthUpdateExpValueForPin(const vector<int> netsID){
    for(size_t pinId = 0;pinId < m_pDB->m_pins.size();pinId++){
        size_t moduleId = m_pDB->m_pins[pinId].moduleId;
        assert(m_usePin[moduleId] == false);// just for test
    }
}

void MyNLP::WireLengthUpdateNetsSumExp(const vector<int> netsID){
    double sum_exp_xi_over_alpha;
    double sum_exp_inv_xi_over_alpha;
    double sum_exp_yi_over_alpha;
    double sum_exp_inv_yi_over_alpha;
    for( unsigned int netId=0; netId < netsID.size(); netId++ ){

        size_t n = netsID[netId];
        if( m_pDB->m_nets[n].size() == 0 )
            continue;
        calc_sum_exp_using_pin(
                    m_pDB->m_nets[n].begin(), m_pDB->m_nets[n].end(), x, _expX,
                    sum_exp_xi_over_alpha, sum_exp_inv_xi_over_alpha,
                    sum_exp_yi_over_alpha, sum_exp_inv_yi_over_alpha );

        m_nets_sum_exp_xi_over_alpha[n]     = sum_exp_xi_over_alpha;
        m_nets_sum_exp_yi_over_alpha[n]     = sum_exp_yi_over_alpha;
        m_nets_sum_exp_inv_xi_over_alpha[n] = sum_exp_inv_xi_over_alpha;
        m_nets_sum_exp_inv_yi_over_alpha[n] = sum_exp_inv_yi_over_alpha;
    }

    for( unsigned int netId=0; netId < netsID.size(); netId++ ){

        size_t n = netsID[netId];
        if( m_pDB->m_nets[n].size() == 0 )
            continue;
        m_nets_sum_p_x_pos[n]     = pow( m_nets_sum_exp_xi_over_alpha[n], 1/_alpha-1 );
        m_nets_sum_p_inv_x_neg[n] = pow( m_nets_sum_exp_inv_xi_over_alpha[n], -1/_alpha-1 );
        m_nets_sum_p_y_pos[n]     = pow( m_nets_sum_exp_yi_over_alpha[n], 1/_alpha-1 );
        m_nets_sum_p_inv_y_neg[n] = pow( m_nets_sum_exp_inv_yi_over_alpha[n], -1/_alpha-1 );

        // not used --hdgao
        /*
        m_nets_sum_p_inv_x_pos[n] = pow( m_nets_sum_exp_inv_xi_over_alpha[n], 1/_alpha-1 );
        m_nets_sum_p_inv_y_pos[n] = pow( m_nets_sum_exp_inv_yi_over_alpha[n], 1/_alpha-1 );
        m_nets_sum_p_x_neg[n]     = pow( m_nets_sum_exp_xi_over_alpha[n], -1/_alpha-1 );
        m_nets_sum_p_y_neg[n]     = pow( m_nets_sum_exp_yi_over_alpha[n], -1/_alpha-1 );
        */

    }
}
// helper
void MyNLP::GetDensityGrad(const size_t moduleId,const vector<int> bin,double &gradx,double &grady){

    double cellX = x[2*moduleId];
    double cellY = x[2*moduleId + 1];
    double width  = m_pDB->m_modules[moduleId].m_width;
    double height = m_pDB->m_modules[moduleId].m_height;
    size_t bin_left,bin_right,bin_bottom,bin_top;
    double left,bottom,right,top;
    if(1){
       double extension_width  = width  * 0.5 + _potentialRX;
       double extension_height = height * 0.5 + _potentialRY;
       left   = cellX - extension_width;
       bottom = cellY - extension_height;
       right  = cellX + extension_width;
       top    = cellY + extension_height;
       if( left   < m_pDB->m_coreRgn.left )	    left   = m_pDB->m_coreRgn.left;
       if( bottom < m_pDB->m_coreRgn.bottom )	bottom = m_pDB->m_coreRgn.bottom;
       if( right  > m_pDB->m_coreRgn.right )	right  = m_pDB->m_coreRgn.right;
       if( top    > m_pDB->m_coreRgn.top )	    top    = m_pDB->m_coreRgn.top;

       int gx, gy;
       GetClosestGrid( left, bottom, gx, gy );

       assert(gx >= 0);
       assert(gy <= m_potentialGridSize);


       bin_left = gx;
       bin_bottom = gy;
       bin_right = m_potentialGridSize - 1;
       bin_top = m_potentialGridSize - 1;
   }
   else{

        bin_left   = bin[0];
        bin_bottom = bin[1];
        bin_right  = bin[2];
        bin_top    = bin[3];

   }
    // test std-cell
    if( height < m_potentialGridHeight && width < m_potentialGridWidth )
        width = height = 0;
    unsigned int gxx,gyy;
    double xx,yy;
    gradx = 0;
    grady = 0;
    for( gxx = bin_left, xx = GetXGrid( bin_left );
           xx < right && gxx <= bin_right;
         gxx++, xx += m_potentialGridWidth ){
        for( gyy = bin_bottom, yy = GetYGrid( bin_bottom );
              yy < top && gyy <= bin_top ;
             gyy++, yy += m_potentialGridHeight ){

            double gX = 0;
            double gY = 0;
            // TEST
            //if( m_gridPotential[ gxx ][ gyy ] > m_expBinPotential[gxx][gyy] )  // TEST for ispd05
            {
                gX = ( m_gridPotential[gxx][gyy] - m_expBinPotential[gxx][gyy] ) *
                        _cellPotentialNorm[moduleId] *
                        GetGradPotential( cellX, xx, _potentialRX, width ) *
                        GetPotential(     cellY, yy, _potentialRY, height );
                gY =  ( m_gridPotential[gxx][gyy] - m_expBinPotential[gxx][gyy] ) *
                        _cellPotentialNorm[moduleId] *
                        GetPotential(     cellX, xx, _potentialRX, width  ) *
                        GetGradPotential( cellY, yy, _potentialRY, height );
            }

            gradx += gX;
            grady += gY;
        }
    } // for each grid
}

void MyNLP::DensityUpdateOneModulePotentialGrid(size_t moduleId,vector<int> bin){
    if( m_pDB->m_modules[moduleId].m_isOutCore )
        return ;
    // preplaced blocks are stored in m_basePotential
    if( m_pDB->m_modules[moduleId].m_isFixed )
        return ;
    // start updating
    double oldcellX = xBak[2*moduleId];
    double oldcellY = xBak[2*moduleId + 1];
    double newcellX = x[2*moduleId];
    double newcellY = x[2*moduleId + 1];
    double width  = m_pDB->m_modules[moduleId].m_width;
    double height = m_pDB->m_modules[moduleId].m_height;
    size_t bin_left,bin_right,bin_bottom,bin_top;
    if(1){
       double left   = newcellX - width  * 0.5 - _potentialRX;
       double bottom = newcellY - height * 0.5 - _potentialRY;
       double right  = newcellX + ( newcellX - left );
       double top    = newcellY + ( newcellY - bottom );
       if( left   < m_pDB->m_coreRgn.left )	left   = m_pDB->m_coreRgn.left;
       if( bottom < m_pDB->m_coreRgn.bottom )	bottom = m_pDB->m_coreRgn.bottom;
       if( right  > m_pDB->m_coreRgn.right )	right  = m_pDB->m_coreRgn.right;
       if( top    > m_pDB->m_coreRgn.top )	    top    = m_pDB->m_coreRgn.top;

       int gx, gy;
       GetClosestGrid( left, bottom, gx, gy );
       assert(gx >= 0);
       assert(gy <= m_potentialGridSize);

       if( gx < 0 )	gx = 0;
       if( gy < 0 )	gy = 0;

       bin_left = gx;
       bin_bottom = gy;
       bin_right = m_potentialGridSize - 1;
       bin_top = m_potentialGridSize - 1;
   }
   else{

    bin_left = bin[0];
    bin_bottom = bin[1];
    bin_right = bin[2];
    bin_top = bin[3];
   }
    // test std-cell
    if( height < m_potentialGridHeight && width < m_potentialGridWidth )
        width = height = 0;
    unsigned int gxx,gyy;
    double xx,yy;
    vector< potentialStruct > potentialList;
    double totalPotential = 0;
    for( gxx = bin_left, xx = GetXGrid( bin_left );
         xx <= bin_right && gxx < m_potentialGridSize;
         gxx++, xx += m_potentialGridWidth ){
        for( gyy = bin_bottom, yy = GetYGrid( bin_bottom );
             yy <= bin_top && gyy < m_potentialGridSize ;
             gyy++, yy += m_potentialGridHeight ){

        }
            double new_potential = GetPotential( newcellX, xx, _potentialRX, width ) *
                    GetPotential( newcellY, yy, _potentialRY, height );
            double old_potential = GetPotential( oldcellX, xx, _potentialRX, width ) *
                    GetPotential( oldcellY, yy, _potentialRY, height );
            potentialList.push_back( potentialStruct( gxx, gyy, new_potential -  old_potential) );
            totalPotential += new_potential;
        }
    double scale = m_pDB->m_modules[moduleId].m_area / totalPotential;
    _cellPotentialNorm[moduleId] = scale;
    for(vector< potentialStruct >::const_iterator ite=potentialList.begin();
        ite!=potentialList.end(); ++ite ) {
        m_gridPotential[ ite->gx ][ ite->gy ] += ite->potential * scale;
    }
}

void MyNLP::GetGrid(int bin_left,int bin_bottom,double bucket_width,double bucket_height,vector<int>& bin){
    double x0 = bin_left*bucket_width + m_pDB->m_coreRgn.left;
    double y0 = bin_bottom*bucket_height + m_pDB->m_coreRgn.bottom;
    double x1 = x0 + bucket_width;
    double y1 = y0 + bucket_height;
    int bin1,bin2,bin3,bin4;
    bin1 = (x0 - m_pDB->m_coreRgn.left)/m_potentialGridWidth;
    bin2 = (y0 - m_pDB->m_coreRgn.bottom)/m_potentialGridHeight;
    bin3 = (x1 - m_pDB->m_coreRgn.left)/m_potentialGridWidth;
    bin4 = (y1 - m_pDB->m_coreRgn.bottom)/m_potentialGridHeight;
    bin3 = (bin3 > m_potentialGridSize - 1) ? m_potentialGridSize - 1 : bin3;
    bin4 = (bin4 > m_potentialGridSize - 1) ? m_potentialGridSize - 1 : bin4;
    bin.push_back(bin1);
    bin.push_back(bin2);
    bin.push_back(bin3);
    bin.push_back(bin4);
}
void MyNLP::GetModulesByNets(const vector<int> netsID,vector<int>& modules){
    vector<bool> isread;
    isread.resize(m_pDB->m_modules.size(),false);
    for(size_t i = 0;i < netsID.size();i++){
        size_t netId = netsID[i];
        size_t moduleSize = m_pDB->m_nets[netId].size();
        for(size_t k = 0;k < moduleSize;k++){
            size_t pinId = m_pDB->m_nets[netId][k];
            size_t moduleId = m_pDB->m_pins[pinId].moduleId;
            if(isread[moduleId]){
                continue;
            }
            else{
                isread[moduleId] = true;
                modules.push_back(moduleId);
            }
        }
    }
}
void MyNLP::WireLengthGradByModule(vector<double>& grad,const vector<int>& modules){
    size_t modulesize = modules.size();
    for(size_t i = 0;i < grad.size();i++)
        grad[i] = 0;
    if( _weightWire > 0 )	//TEST
        for( unsigned int k = 0; k < modulesize; k++ )	// for each block
        {
            size_t i = modules[k];
            if( m_pDB->m_modules[i].m_isFixed || m_pDB->m_modules[i].m_netsId.size() == 0 )
                continue;



            for( unsigned int j=0; j<m_pDB->m_modules[i].m_netsId.size(); j++ )
            {
                // for each net connecting to the block
                size_t netId = m_pDB->m_modules[i].m_netsId[j];
                if( m_pDB->m_nets[netId].size() == 0 ) // floating-module
                    continue;

                // TODO: modification for LEF/DEF input
                // no floating pin for bookshelf format
                //if( m_pDB->m_nets[netId].size() == 0 )
                //	continue;

                int selfPinId = m_moduleNetPinId[i][j];

                if( m_usePin[i] )
                {
                    assert( selfPinId != -1 );
                    double xx = x[ 2*i ]   + m_pDB->m_pins[ selfPinId ].xOff;
                    double yy = x[ 2*i+1 ] + m_pDB->m_pins[ selfPinId ].yOff;
                    xx *= m_posScale;
                    yy *= m_posScale;

                    grad[ 2*i ] +=
                            m_nets_sum_p_x_pos[netId]     * _expPins[2*selfPinId] / xx -
                            m_nets_sum_p_inv_x_neg[netId] / _expPins[2*selfPinId] / xx;
                    grad[ 2*i+1 ] +=
                            m_nets_sum_p_y_pos[netId]     * _expPins[2*selfPinId+1] / yy -
                            m_nets_sum_p_inv_y_neg[netId] / _expPins[2*selfPinId+1] / yy;
                }
                else
                {
                    double xx = x[ 2*i ];
                    double yy = x[ 2*i+1 ];
                    xx *= m_posScale;
                    yy *= m_posScale;

                    grad[ 2*i ] +=
                            m_nets_sum_p_x_pos[netId]     * _expX[2*i] / xx  -
                            m_nets_sum_p_inv_x_neg[netId] / _expX[2*i] / xx;
                    grad[ 2*i+1 ] +=
                            m_nets_sum_p_y_pos[netId]     * _expX[2*i+1] / yy -
                            m_nets_sum_p_inv_y_neg[netId] / _expX[2*i+1] / yy;
                }

            } // for each pin in the module
        } // for each module
    for(size_t i = 0;i < grad.size();i++)
            grad[i] *= _weightWire;
}
void MyNLP::DensityUpdatePotentialGridByModule(const vector<int>& modules,const vector<int>& bin ){
    size_t moduleSize = modules.size();
    for(size_t i = 0;i < moduleSize;i++){
        size_t moduleId = modules[i];
        UpdateModulePotentialGridByDiff(moduleId,bin);
    }
/*
    for(size_t i = 0;i < m_potentialGridSize;i++)
        for(size_t j = 0;j < m_potentialGridSize;j++){
            if( m_gridPotential[i][j] < 0 )
                m_gridPotential[i][j] = 0;
        }
*/
}
void MyNLP::UpdateModulePotentialGridByDiff(size_t moduleId,vector<int> bin){
    if( m_pDB->m_modules[moduleId].m_isOutCore )
        return ;
    // preplaced blocks are stored in m_basePotential
    if( m_pDB->m_modules[moduleId].m_isFixed )
        return ;
    double oldcellX = xBak[2*moduleId];
    double oldcellY = xBak[2*moduleId + 1];
    double newcellX = x[2*moduleId];
    double newcellY = x[2*moduleId + 1];
    double width  = m_pDB->m_modules[moduleId].m_width;
    double height = m_pDB->m_modules[moduleId].m_height;
    size_t bin_left,bin_right,bin_bottom,bin_top;
    double left,right,bottom,top;
    if(1){

        double extension_width  = width*0.5  + _potentialRX;
        double extension_height = height*0.5 + _potentialRY;

        if(oldcellX > newcellX){
            left  = newcellX - extension_width;
            right = oldcellX + extension_width;
        }
        else{
            left  = oldcellX - extension_width;
            right = newcellX + extension_width;
        }
        if(oldcellY > newcellY){
            bottom = newcellY - extension_height;
            top    = oldcellY + extension_height;
        }
        else{
            bottom = oldcellY - extension_height;
            top    = newcellY + extension_height;
        }
#if 0
        left   = newcellX - extension_width;
        right  = newcellX + extension_width;
        bottom = newcellY - extension_height;
        top    = newcellY + extension_height;
#endif
        if( left   < m_pDB->m_coreRgn.left )	left   = m_pDB->m_coreRgn.left;
        if( bottom < m_pDB->m_coreRgn.bottom )	bottom = m_pDB->m_coreRgn.bottom;
        if( right  > m_pDB->m_coreRgn.right )	right  = m_pDB->m_coreRgn.right;
        if( top    > m_pDB->m_coreRgn.top )	    top    = m_pDB->m_coreRgn.top;

        int gx,gy;
        GetClosestGrid( left, bottom, gx, gy );
        assert(gx >= 0);
        assert(gy <= m_potentialGridSize);

        bin_left   = gx;
        bin_bottom = gy;
        bin_right  = m_potentialGridSize - 1;//not used
        bin_top    = m_potentialGridSize - 1;//not used
    }
    else{
        bin_left   = bin[0];
        bin_bottom = bin[1];
        bin_right  = bin[2];
        bin_top    = bin[3];
    }
#if DEBUG_hd
    printf("moduleId = %d,x1 = %f,x2 = %f,left = %f,bottom = %f\n",moduleId,newcellX,newcellY,left,bottom);
#endif
    if( height < m_potentialGridHeight && width < m_potentialGridWidth )
        width = height = 0;
    unsigned int gxx,gyy;
    double xx,yy;
    vector< potentialStruct > potentialList;
    double totalPotential = 0;
    for( gxx = bin_left, xx = GetXGrid( bin_left );xx <= right && gxx < bin_right;gxx++, xx += m_potentialGridWidth ){
        for( gyy = bin_bottom, yy = GetYGrid( bin_bottom );yy <= top && gyy < bin_top ;gyy++, yy += m_potentialGridHeight){
            double new_potential = GetPotential( newcellX, xx, _potentialRX, width ) *
                    GetPotential( newcellY, yy, _potentialRY, height );
            double old_potential = GetPotential( oldcellX, xx, _potentialRX, width ) *
                    GetPotential( oldcellY, yy, _potentialRY, height );
            if(new_potential < 0) new_potential = 0;
            if(old_potential < 0) old_potential = 0;
            double delta = new_potential -  old_potential;
            potentialList.push_back( potentialStruct( gxx, gyy, delta) );
            totalPotential += new_potential;
#if DEBUG_hd
            printf("gx = %-4d,gy = %-4dnew  = %-4.3f,old = %-4.3f,delta = %-4.3f\n",gxx,gyy,new_potential,old_potential,delta);
            fflush(stdout);
            assert(new_potential == new_potential);
            assert(old_potential == old_potential);
            assert(delta == delta);
#endif
        }
    }
#if DEBUG_hd
    printf("totalPotential = %-4.3f\n",totalPotential);
    assert(totalPotential != 0.0);
#endif
    double scale = m_pDB->m_modules[moduleId].m_area / totalPotential;
    _cellPotentialNorm[moduleId] = scale;// minor bug here
    for(vector< potentialStruct >::const_iterator ite=potentialList.begin();
        ite!=potentialList.end(); ++ite ) {
        m_gridPotential[ ite->gx ][ ite->gy ] += (ite->potential * scale);
    }

}
/********************************************hdgao added***************************************************************/

