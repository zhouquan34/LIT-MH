#include "model.h"

using namespace std;

int g_init_add_mark = -1;
int g_init_add_count = 0; 
vector<double> g_proposal_weights; 
vector<int> g_pass; 

double sqd (double x, double y){
	return (x-y)*(x-y);
}

model::model(void){
	sse = 0.0; pi = 0.0;  pve_rb = -1.0;
	tau = UNDEF_LL; like = 0.0; 
	gamma_rb_mse = 0.0; beta_rb_mse = 0.0;		
	mode = -1; 
	
	R = NULL;  Beta = NULL; 
	XX = NULL; X = NULL;    Xy = NULL;
	Rinv = NULL; Rinv_xy = NULL;
	
//	if (g_rw == 0){
		neighbor_ll.assign(g_n_snp, UNDEF_LL);
		neighbor_beta.assign(g_n_snp, UNDEF_LL);
//	}else{
//		neighbor_ll.assign(g_n_snp, 0);
//		neighbor_beta.assign(g_n_snp, 0);		
//	}
}

void model::clean(void){
	free1D(R);	
	free1D(XX);	
	free1D(Xy);
	free1D(Beta);
	free1D_pointer(X);
	free1D(Rinv);
	free1D(Rinv_xy);
}

model::~model(){
	clean();
}

double model::model_ll(void){return like;}

int model::model_size(void){return subset.size();}

int model::find_snp(int d){
	for (int j=0; j<subset.size(); j++){
		if (subset[j] == d){ return j; }
	}
	return -1; 
}

string model::model_str(void){
	stringstream s;
	if (subset.size() == 0){ return s.str(); }
	s << subset[0];
	for (int i=1; i<subset.size(); i++){ s << "," << subset[i]; }
	return s.str();
}

string model::para_str(void){
	stringstream s;
// assuming no covariate	
	double* hatY = allocate1D(g_n_sub);
	double* z = allocate1D(msize);
	for (int i=0; i<msize; i++){z[i] = gsl_ran_gaussian(gsl_r, 1.0/sqrt(tau));}
	double* beta_sampled = allocate1D(msize);
	Chol_solve_backward(msize, R, z, beta_sampled);
	for (int i=0; i<msize; i++){beta_sampled[i] += Beta[i + g_n_cov];}
	mat_vec_mul(X, beta_sampled, msize, g_n_sub, hatY, 1);
	double pve = array_sst(hatY, g_n_sub) / g_yy;
	free1D(z); 
	free1D(hatY); 
	free1D(beta_sampled);

	s << subset.size() << "\t"  << tau << "\t" << like << "\t" << pve << "\t" << pve_rb;
	if (g_true == 1){
		double e1 = 0;
		double e2 = 0;
		vector<double> bb (g_n_snp, 0);
		for (int k=0; k<subset.size(); k++){bb[subset[k]] = Beta[k + g_n_cov];}
		for (int i=0; i<g_n_snp; i++){
			e1 += sqd(selected[i], g_true_gamma[i]);
			e2 += sqd(bb[i], g_true_beta[i]);
		}	
		s << "\t" << e1 << "\t" << gamma_rb_mse << "\t" << e2 << "\t" << beta_rb_mse;		
	}
	
	if (mode > -1){s << "\t" << mode;}
	return s.str();
}


void model::update_marginal(int k){
	for (int i=0; i<subset.size(); i++){
		g_pip[subset[i]] += k;
		g_beta[subset[i]] += k * Beta[i + g_n_cov];
	}
	g_model_size[subset.size()] += k;
}

void model::allocate(int p){
	clean();
	X = allocate1D_float_pointer(p);
	XX = allocate1D(p * (p+1) / 2);
	R = allocate1D(p * p);
	Xy = allocate1D(p);
	Beta = allocate1D(p);
}

void model::initialize(vector<int>& start){
	subset = start;
	msize = subset.size() + g_n_cov;
	allocate(msize);	
	X_init(subset, X);
	XX_init(subset, X, msize, XX);
	Xy_init(X, msize, Xy);
	Chol_init(msize, XX, R);
	
	selected.clear();
	selected.assign(g_n_snp,0);
	for (int i=0; i<subset.size(); i++){ selected[subset[i]] = 1; }
	
	calc_rss();
}

void model::copy(model* m){  
	pi = m->pi;
	sse = m->sse;
	like = m->like;
//	tau = m->tau; 
	msize = m->msize;
	subset = m->subset;
	selected = m->selected;
}


void model::init_add(model* m, vector<int> add){	
	int add_size = add.size(); 
	copy(m); 
	msize = m->msize + add_size;
	allocate(msize); 

	X_add(add, m->msize, m->X, X); 		
	XX_add(subset, add, m->msize, X, m->XX, XX); 
	Xy_add(add, m->msize, m->Xy, Xy);
	Chol_add(1, m->msize, XX, m->R, R); 
	
	for (int i=0; i<add.size(); i++){
		selected[add[i]] = 1; 
		subset.push_back(add[i]);
	}
	
	calc_rss(); 
}

void model::init_add(model* m, int k){
	vector<int> add(1, k);
	init_add(m, add);
}

void model::init_del(model* m, vector<int> del){
	sort(del.begin(), del.end()); 
	int del_size = del.size();
	copy(m);
	msize = m->msize - del_size; 
	allocate(msize);

	X_del(del,  m->msize, m->X,  X); 		
	XX_del(del, m->msize, m->XX, XX); 
	Xy_del(del, m->msize, m->Xy, Xy);
	Chol_del(del, m->msize, m->R, R);

	for (int i=del.size() - 1; i>=0; i--){
		selected[subset[del[i]]] = 0; 
		subset.erase( subset.begin() + del[i] );
	}
	
	calc_rss(); 
}	

void model::init_del(model* m, int k){
	vector<int> del(1, k);
	init_del(m, del);
}

void model::sample_tau(void){
	if (tau != UNDEF_LL){return;}
	tau = gsl_ran_gamma(gsl_r, g_n_sub / 2.0, 2.0 / sse);	
}

// like actually means unnormalized log-posterior. model size penalty is considered. 
void model::calc_rss(void){	
	Chol_solve(msize, R, Xy, Beta);		
	sse = g_yy - vec_vec_mul(Xy, Beta, msize); 
	like = -0.5 * g_n_sub  * log(sse / g_yy + g_gg_const) - subset.size() * 0.5 * g_gg_exp * log(g_n_snp);
	
	if (g_use_kappa == 1){
		like = like - subset.size() * g_kappa * log(g_n_snp) ;
	}else{
		like = like - subset.size() * log(pi);
	}
	sample_tau();
}

double model::add_ll(int k){
	double q = 0;
	if (neighbor_ll[k] != UNDEF_LL){
		q = neighbor_ll[k];
	}else{
	//	model* nb_model = new model();
	//	nb_model->init_add(this, k); 		
	//	q = nb_model->like - like;
	//	neighbor_ll[k] = q;
	//	neighbor_beta[k] = nb_model->Beta[g_n_cov + subset.size()];
	//	delete(nb_model); 
	
		// assume g_n_cov = 0;
		if (!Rinv){cerr << "Wrong: Rinv not calculated yet" << endl;}
		double* Rinv_xj = allocate1D(msize);
		double* Xj = allocate1D(msize); 
		for (int i=0; i<msize; i++){
			Xj[i] = XX_get(k, subset[i]);
		}
		for (int i=0; i<msize; i++){
			double s = 0;
			for (int j=0; j<=i; j++){
				s += Rinv[j * msize + i] * Xj[j];
			}
			Rinv_xj[i] = s; 
		}
		double t1 = g_xty[k] - vec_vec_mul(Rinv_xj, Rinv_xy, msize);
		double t2 = g_snp_var[k] * g_n_sub - vec_vec_mul(Rinv_xj, Rinv_xj, msize);
		q = -0.5*g_n_sub*log(1.0 - t1*t1/t2/(sse + g_yy*g_gg_const)) - 0.5*g_gg_exp*log(g_n_snp);
		if (g_use_kappa == 1){q = q - g_kappa * log(g_n_snp);}
		else{q = q - log(pi);}
		neighbor_ll[k] = q;
		//cout << k << " " << q0 << " " << q << " " << sse0 << " " << sse - t1*t1/t2  << " " << sse << " " << t1*t1/t2 << endl;
		neighbor_beta[k] = t1/t2; 
		gsl_ran_gamma(gsl_r, g_n_sub / 2.0, 2.0 / sse);	// reproduce the previous results calcuated using chol.  
		free1D(Xj);
		free1D(Rinv_xj);			
	}
	return q; 	
}

double model::add_weight(int k){
	if (g_rw == 1){return 0;}
	if (g_pass.size() > 0){
		if (g_pass[k] == 1){return g_lower_add * log((double) g_n_snp);}
	}		

	double q = 0;
	if (g_approx_informed == 1){q = g_proposal_weights[k];}
	else { q  = add_ll(k); }

	if (g_lbsq == 1){return q/2;}

	double lp = log((double) g_n_snp);
	if (q > g_upper_add * lp ){ q = g_upper_add * lp ;}
	if (q < g_lower_add * lp ){ q = g_lower_add * lp ;}
	return q; 
}

double model::del_ll(int k){
//	if (g_rw == 1){return 0;}
	
	int s = subset[k];	
	double q = 0;
	if (neighbor_ll[s] != UNDEF_LL){
		q = neighbor_ll[s];
	}else{
		model* nb_model = new model();
		nb_model->init_del(this, k); 		
		q = nb_model->like - like;
		//cout << "D2: " << nb_model->like << "\t" << like << endl;
		neighbor_ll[s] = q;
		delete(nb_model); 
	}

	return q; 	
}


double model::del_weight(int k){
	if (g_rw == 1){return 0;}
	double q = 0;
	if (g_approx_informed == 1){q = -g_proposal_weights[subset[k]];}
	else {q = del_ll(k); }	

	if (g_lbsq == 1){return q/2;}
		
	double lp = log((double) g_n_snp);
	if (q > g_upper_del * lp ){ q = g_upper_del * lp;}
	if (q < g_lower_del * lp ){ q = g_lower_del * lp;}
	return q; 
}

void model::precalc_add(){
	if (Rinv){return;}	

	record_time();	
	Rinv = allocate1D(msize * msize);		
	for (int j = 0; j < msize; j++){
		Rinv[j*msize + j] = 1.0 / R[j*msize + j];
		for (int i = j - 1; i >= 0; i--){
			double s = 0.0;
			for (int k = j ; k > i; k --){
				s += Rinv[k * msize + j] * R[i * msize + k];
			}
			Rinv[i*msize + j] =  -s / R[i*msize + i];
		}
	}

	Rinv_xy = allocate1D(msize);
	for (int i=0; i<msize; i++){
		double s = 0;
		for (int j=0; j<=i; j++){
			s += Rinv[j * msize + i] * Xy[j];
		}
		Rinv_xy[i] = s; 
	}
	record_time("Precalc_add");
//	cout << "R: " << msize << endl;
//	print_matrix(R, msize*msize, msize);
//	cout << "Rinv: " << msize << endl;
//	print_matrix(Rinv, msize*msize, msize);
	return;
}


// if_swap = 1: no effect on adding. 
int model::propose_add(double& log_prop){
	if (g_rw == 0){ precalc_add();}
	
	vector<int> ids; 
	vector<double> pp; 	
	for (int i=0; i<g_n_snp; i++){
		if (selected[i] == 0){
			ids.push_back(i);
			pp.push_back(add_weight(i));
		}
	}
	int k = log_weighted_sample(pp, log_prop);
	return ids[k];
}

double model::calc_add(int k){
	if (g_rw == 0){ precalc_add();}
	
	double s = -1e10;
	double qk = 0;
	int found = 0; 
	for (int i=0; i<g_n_snp; i++){
		if (selected[i] == 0){
			double q = add_weight(i);
			s = log_sum_log(s, q);
			if (i == k){ qk = q; found = 1;}
		}
	}
	if (found == 0){cerr << "Warning! variable to add does not exist!" << endl;}
	return qk - s; 
}
	
// if swap, cannot remove the last one. 
int model::propose_del(double& log_prop, int if_swap){
	vector<double> pp; 	
	for (int j=0; j<subset.size() - if_swap; j++){
		pp.push_back(del_weight(j));
	}
	//cout << "del" << endl; print_vector(pp); cout << endl; 
	return log_weighted_sample(pp, log_prop);
}

// avoid default -1
double model::calc_del(int k, int avoid){
	//print_vector(subset);
	//cout << "calc_del: " << k << endl;
	double s = -1e10;
	double qk = 0;
	int found = 0;	
	for (int j=0; j<subset.size(); j++){
		if (j == avoid){continue;}
		double q = del_weight(j);
		s = log_sum_log(s, q); 
		//cout << "D1: " << j << "\t" << q << "\t" << s  << "\t" << k << endl;
		if (j == k){qk = q; found = 1;}
	}
	
	if (found == 0){cerr << "Warning! variable to delete does not exist! " << k << "\t" << avoid << endl;}
	return qk - s;
}

double model::add_snp (model* last_model){	
	if (last_model->subset.size() + 1 > g_max_model_size){return FAIL_ALPHA;} 
	double log_prop_new = 0.0;
	int add = last_model->propose_add(log_prop_new); 
	init_add(last_model, add);	
	double log_prop_old = calc_del(subset.size() - 1);
	return log_prop_old - log_prop_new;
}

double model::del_snp (model* last_model){
	if (last_model->subset.size() < (1 + g_min_model_size) ){return FAIL_ALPHA;}	
	double log_prop_new = 0.0; 
	int del = last_model->propose_del(log_prop_new);	
	init_del(last_model, del); 
	double log_prop_old = calc_add(last_model->subset[del]); 
	return log_prop_old - log_prop_new;
}

int model::if_swap_allowed(void){
	if (subset.size() + 1 > g_max_model_size){return 0;} 
	if (subset.size() < (1 + g_min_model_size) ){return 0;}	
	return 1;
}

double model::swap_snp (model* last_model){
	if (last_model->subset.size() + 1 > g_max_model_size){return FAIL_ALPHA;} 
	if (last_model->subset.size() < (1 + g_min_model_size) ){return FAIL_ALPHA;}	

	double log_prop_new_add = 0.0;
	double log_prop_new_del = 0.0;
	model* tmp_model = new model();
	int add = last_model->propose_add(log_prop_new_add); 
	tmp_model->init_add(last_model, add);	
	int del = tmp_model->propose_del(log_prop_new_del, 1);	
	init_del(tmp_model, del); 
	
	double log_prop_old_del = tmp_model->calc_del(tmp_model->subset.size() - 1, del);
	double log_prop_old_add = calc_add(tmp_model->subset[del]); 
	delete(tmp_model);
	
//	cout << "swap: " << last_model->subset.size() << "\t" << log_prop_old_add << "\t" << log_prop_old_del << "\t" << log_prop_new_add << "\t" << log_prop_new_del << endl;
	
	return log_prop_old_add + log_prop_old_del - log_prop_new_add - log_prop_new_del;
}


int model::compare_models(double alpha, model* last_model){
	double mh = alpha + (like - last_model->like);
	double r = gsl_rng_uniform(gsl_r);

	if ( mh > log(r)){
		return 1; // accept
	}else{
		return 0; // reject
	}
}

// needs to be modified if g_n_cov > 0
void model::rao_blackwell (void){
	precalc_add(); // should be prohibited when using rw
	record_time();
	
	sample_tau();	
	double* beta_rb = allocate1D(g_n_snp);
	//copy1D(beta_rb, Beta, g_n_cov);
	double correct = 0.0; 	// approximate  
	gamma_rb_mse = 0.0;	
	beta_rb_mse = 0.0;	


	for (int i=0; i<g_n_snp; i++){
		if (selected[i] == 0){
			double wt = add_ll(i);
			g_proposal_weights[i] = wt; 
			double q = neighbor_ll[i];
			double u = log_sum_log(0, q);
			double pr = exp(q - u);
			double b0 = neighbor_beta[i];
			g_rb_pip[i] += pr;
			gamma_rb_mse += sqd(pr, g_true_gamma[i]);
			double b = b0 * pr;
			g_rb_beta[i] += b;
			beta_rb[i] = b;
			beta_rb_mse += sqd(b, g_true_beta[i]);
			double vxx = 1.0/tau/g_snp_var[i]; 
			double var_beta_rb = pr*(b0*b0 + vxx) - pr*pr*b0*b0;
			correct += g_snp_var[i] * var_beta_rb;
		}
	}


	for (int j=0; j<subset.size(); j++){
		int k = subset[j];
		double wt = del_ll(j);
		g_proposal_weights[k] = -wt; 
		double u = log_sum_log(0, neighbor_ll[k]);
		double pr = exp(-u);
		g_rb_pip[k] += pr;
		gamma_rb_mse += sqd(pr, g_true_gamma[k]);
		double b0 = Beta[g_n_cov + j];
		double b = b0 * pr;
		g_rb_beta[k] += b;
		//beta_rb[k + g_n_cov] = b;
		beta_rb[k] = b;
		beta_rb_mse += sqd(b, g_true_beta[k]);
		double vxx = 1.0/tau/g_snp_var[k]; 
		double var_beta_rb = pr*(b0*b0 + vxx) - pr*pr*b0*b0;
		correct += g_snp_var[k] * var_beta_rb;
	}

	
	double* yhat = allocate1D(g_n_sub);
	double* resid = allocate1D(g_n_sub);
	mat_vec_mul(g_data_mat, beta_rb, g_n_snp, g_n_sub, yhat, 1);
	double ssr = array_sst(yhat, g_n_sub);
	pve_rb = ssr / g_yy; 
	for (int i=0; i<g_n_sub; i++){resid[i] = g_pheno[i] - yhat[i];}
	double sse = array_sst(resid, g_n_sub);
	free1D(yhat);	
	free1D(resid); 
	pve_rb = (ssr + correct) / g_yy; // this is also how we calculcate the pve in simulation. 
	free1D(beta_rb);

	mode = 1;
	for (int i=0; i<g_n_snp; i++){
		if (neighbor_ll[i] > 0){mode = 0; break;}
	}

	g_rb_count ++;
	record_time("RaoBlackwell");
	return;
} 

