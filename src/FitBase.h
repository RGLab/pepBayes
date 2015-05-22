
#ifndef FITBASE_H_
#define FITBASE_H_

static inline void check_interrupt_impl(void* /*dummy*/) {
    R_CheckUserInterrupt();
}

inline bool check_interrupt() {
    return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}

class EcmAlgorithm {
public:
    int n_iter, n_threads;

    EcmAlgorithm() :
        n_iter(10),
        n_threads(1) {};

    EcmAlgorithm(const Rcpp::List & iter_list) {
        n_iter = Rcpp::as<int>(iter_list["n_iter"]);
        n_threads = Rcpp::as<int>(iter_list["n_threads"]);
    };

    virtual ~EcmAlgorithm() {};

    Rcpp::List run() {
        bool input_ok = checkInput();
        bool interrupt = false;
        if (!input_ok)
            Rcpp::stop("invalid input");
        try {
            for(int i = 0; (i < n_iter) & (!interrupt); i++) {
                iterate(i);
                collectIteration(i);
                if (check_interrupt())
                    interrupt = true;
            }
            return gatherOutput();
        } catch (...) {
            return gatherOutput();
            //::Rf_error("c++ exception (unknown reason)");
        }
        Rcpp::List out = Rcpp::List::create(Rcpp::Named("") = R_NilValue);
        return out;
    }

private:
    virtual bool checkInput() =0;
    virtual void iterate(const int & iter_idx) =0;
    virtual void collectIteration(const int & iter_idx) =0;
    virtual const Rcpp::List gatherOutput() =0;
};

class MHTuner {
    int num_accept;
    int num_trials;
    int total_trials;
    int n_burn;
    int adjustment_interval;
    double target_prop;
    double scale;
    void update_scale() {
        const double emp_prop = ((double) num_accept) / num_trials;
        if (emp_prop < target_prop - .05) {
            scale *= 1.0 - .4 * (target_prop - emp_prop) / target_prop;
        } else if (emp_prop > target_prop + .05) {
            scale *= 1.0 + .4 * (emp_prop - target_prop) / (1.0 - target_prop);
        }
        num_accept = 0;
        num_trials = 0;
    }
public:
    MHTuner() : num_accept(0), num_trials(0), total_trials(0), n_burn(0),
            adjustment_interval(100), target_prop(.40),
            scale(log(2)) {};

    MHTuner(int burn) : num_accept(0), num_trials(0), total_trials(0), n_burn(burn),
            adjustment_interval(100), target_prop(.40),
            scale(log(2)) {};

    void update(bool accepted) {
        num_accept += accepted;
        num_trials++;
        total_trials++;
        if (num_trials == adjustment_interval && num_trials < n_burn)
            update_scale();
    };

    double getScale() const {
        return scale;
    }
};

class MarkovChain {
public:
    int n_samples, n_burn, n_thin, n_threads, seed;
    MarkovChain(const Rcpp::List & chain_list) :
        n_samples(Rcpp::as<int>(chain_list["n_samples"])),
        n_burn(Rcpp::as<int>(chain_list["n_burn"])),
        n_thin(Rcpp::as<int>(chain_list["n_thin"])),
        n_threads(Rcpp::as<int>(chain_list["n_threads"])),
        seed(Rcpp::as<int>(chain_list["seed"])) {};

    virtual ~MarkovChain() {};

    Rcpp::List run() {
        // burn in the chain
        bool interrupt = false;
        for(int i = 0; (i < n_burn) & (!interrupt); i++) {
            // chains should NOT implement both of these
            // or else...
            iterate();
            if (check_interrupt())
                interrupt = true;
        }
        // collect samples
        int sample_idx = 0;
        const int total_iter = n_samples * n_thin + n_burn;
        for(int i = n_burn + 1; (i <= total_iter) & (!interrupt); i++) {
            iterate();
            if ((i - n_burn) % n_thin == 0) {
                if (check_interrupt())
                    interrupt = true;
                collectIteration(sample_idx);
                sample_idx++;
            }
        }
        if (interrupt)
            return Rcpp::List::create(Rcpp::Named("NULL") = R_NilValue);
        return chainOutput();
    }
private:
    virtual void iterate() {};
    virtual void collectIteration(const int & sample_idx) =0;
    virtual const Rcpp::List chainOutput() =0;
};

#endif
