#include "mcmc.h"

using namespace std;

void proposal_type (int& add, int& del, int& prop_type, int if_swap){	
	if (gsl_rng_uniform(gsl_r) < g_long_range){ 
		prop_type = 2;
		int csize = gsl_rng_uniform_int(gsl_r, g_max_jump) + 1; 
		add = csize; del = csize;
		//cout << "long-range " << csize << endl; 
	}else{
		if (gsl_rng_uniform(gsl_r) > g_prop_add){
			prop_type = 1;
			add = 0; del = 1;
		}else{
			prop_type = 0;
			add = 1; del = 0; 
		}
	}
}

void mcmc_sample(model*& last_model){
	int add = 0, del = 0;
	int prop_type = 0;
	//proposal_type(add, del, prop_type, last_model->if_swap_allowed() );
	proposal_type(add, del, prop_type, 1 );
	//cout << g_mcmc_i << " type: " << prop_type << "; add: " << add << "; del:" << del << endl;

	double alpha = 0.0;
	model* proposed_model = new model();
	int if_last = 0;
	if (prop_type == 0){ 
		alpha = proposed_model->add_snp(last_model);	
	}else if (prop_type == 1){  
		alpha = proposed_model->del_snp(last_model);
	}else{ // long-range
		alpha = proposed_model->swap_snp(last_model);
	/*	model* tmp_model = new model();
		for (int k=1; k<add; k++){
			alpha += proposed_model->swap_snp(tmp_model);
		}
		delete(tmp_model);
	*/
	}

	//cout << "alpha: " << alpha << endl;

	int acc = 0; 
	if (alpha != FAIL_ALPHA){
		acc = proposed_model->compare_models(alpha, last_model);
	}

	if (acc == 1){	
		if (g_mcmc_i > 0 ){
			string paras = last_model->para_str();
			string model = last_model->model_str();
			PATH << g_mcmc_i - 1 << "\t" << g_mcmc_stay << "\t"  << paras << endl;
			MOD << model << endl;
		}
		if (g_mcmc_i > g_mcmc_warm) {
			int s = g_mcmc_i - g_mcmc_warm - 1; 
			if (s > g_mcmc_stay) {s = g_mcmc_stay;}
			last_model -> update_marginal(s);
		}
		delete(last_model);
		last_model = proposed_model;
		g_mcmc_stay = 1;
	}else{
		delete(proposed_model);
		g_mcmc_stay ++ ;
	}
		
	if (g_mcmc_i > g_mcmc_warm){
		g_propose[prop_type] ++;
		if (acc == 1) {g_accept[prop_type] ++;}
	}	

	return;
}


void mcmc_init(model*& model){
	if (g_max_model_size > g_n_snp) {g_max_model_size = g_n_snp;}
	if (g_min_model_size < 0){g_min_model_size = 0;}
	if ((g_upper_add - g_lower_add < 1e-3) && (g_upper_del - g_lower_del < 1e-3)){
		g_rw = 1; 
	}
	// this is commented for counting the number of local modes
	if (g_rw == 1){
		g_rb_thin = g_mcmc_iter + g_mcmc_warm + g_last_iter + 100;
		cerr << "Use random walk proposals and no Rao-Blackwellization" << endl;	
	}

	if (g_skip == 1){
		g_rb_thin = g_mcmc_iter + g_mcmc_warm + g_last_iter + 100;
		cerr << "No Rao-Blackwellization" << endl;	
	}
	//if (g_exact_bf == 1) {g_icf_min_p = g_n_snp + 1;}
	
	g_pip.clear(); g_pip.assign(g_n_snp, 0);
	g_rb_pip.clear(); g_rb_pip.assign(g_n_snp, 0.0);
	g_beta.clear(); g_beta.assign(g_n_snp, 0.0);
	g_rb_beta.clear(); g_rb_beta.assign(g_n_snp, 0.0);	
	g_propose.assign(3,0);
	g_accept.assign(3,0);

	if (g_init_model.size() > 0){  
		model->initialize(g_init_model);
		LOG << "Starting model size = " << g_init_model.size() << endl;
	}else{
		vector<int> s;	
		if (g_start_size >= 0){
			vector<int> ids;
			for (int i=0; i<g_n_snp; i++){ ids.push_back(i); }
			if (g_ordered_start == 0){
				random_shuffle(ids.begin(), ids.end());
			}
			for (int i=0; i<g_start_size; i++){
			//	s.push_back(g_single_id[i]);
				s.push_back(ids[i]);
			}
		}else{
			for (int i=0; i<g_n_snp; i++){
				if (g_true_gamma[i] == 1){
					s.push_back(i);
				}
			}	
		}
		model->initialize(s);
		LOG << "Starting model size = " << s.size() << endl;
	}
	
	if (g_rw == 0 && g_skip == 0){model->rao_blackwell();}	
	g_chol_call=0;

	return;
}

int calc_refresh(void){
	int u = 0;
	if (g_continue == 0){ u = (g_mcmc_warm + g_mcmc_iter) / 50; }
	else{ u = g_mcmc_iter / 50; }
	if (u <= 0){u = 1;}
	return u;
}

void calc_start_end(int& s, int& e){
	if (g_continue == 0){
		s = 1;
		e = g_mcmc_iter + g_mcmc_warm;
	}else{
		s = g_last_iter + g_mcmc_warm + 1;
		e = g_mcmc_iter + g_mcmc_warm + g_last_iter;
	}
}

void mcmc (void){
	cerr << "Initializing MCMC" << endl;
	g_diag_c = 0.0; 
	g_proposal_weights.assign(g_n_snp, 0);
	model* my_model = new model(); 
	mcmc_init(my_model);
	
	int print_unit = calc_refresh();
	g_mcmc_stay = 0;
	cerr << "MCMC sampling" << endl;
	if (g_true == 0){
		PATH << "No.\tLife\tModelSize\ttau\tlog-like\tpve\tpve-rb" << endl;
	}else{
		PATH << "No.\tLife\tModelSize\ttau\tlog-like\tpve\tpve-rb\tgamma-se\tgamma-rb-se\tbeta-se\tbeta-rb-se" << endl;
	}
	MOD << "#Iteration index matched with the PATH file" << endl;
	
	int start_iter = 0, total_iter = 0;
	calc_start_end(start_iter, total_iter);
	for (g_mcmc_i = start_iter; g_mcmc_i <= total_iter; g_mcmc_i ++){
		if (g_mcmc_i % print_unit == 0){cerr << '='; fflush(stdout);}
		mcmc_sample(my_model);
		if (g_mcmc_i > g_mcmc_warm){
			if (g_mcmc_i % g_rb_thin == 0){ my_model->rao_blackwell();}
		}
	}
	
	string paras = my_model->para_str();
	string model = my_model->model_str();
	PATH << g_mcmc_i - 1  << "\t" << g_mcmc_stay  << "\t" << paras << endl;	
	MOD << model << endl;
	my_model -> update_marginal(g_mcmc_stay);
	
	delete(my_model);	
	mcmc_output();	
	
	return;
}	


