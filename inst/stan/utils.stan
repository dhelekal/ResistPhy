    real interp_seg_lin(real t, real t_lower, real t_upper, real usage_lower, real usage_upper) {
        real slope = (usage_upper - usage_lower)/(t_upper - t_lower);
        real s = t - t_lower;
        real out = s*slope + usage_lower;
        return out;
    }  

    vector interp_usage_step(vector times, array[] real ABX_usage, array[] real usage_ts, int N, int M) {
        vector[N] out;
        for (i in 1:N) {
            real t = times[i];
            for(j in 1:(M-1)) {
                if (t >= usage_ts[j] && t < usage_ts[j+1]){
                    out[i] = ABX_usage[j];
                } else if (t >= usage_ts[M]) {
                    out[i] = ABX_usage[M];
                }
            }
        }
        return out;
    } 

    vector interp_usage_lin(vector times, array[] real ABX_usage, array[] real usage_ts, int N, int M) {
        vector[N] out;
        for (i in 1:N) {
            real t = times[i];
            if(t == usage_ts[1]) {
                out[i] = ABX_usage[1];
                continue;
            }
            for(j in 2:M) {
                if (t > usage_ts[j-1] && t <= usage_ts[j]){
                    out[i] = interp_seg_lin(t, usage_ts[j-1], usage_ts[j], ABX_usage[j-1], ABX_usage[j]);
                }
            }
        }
        return out;
    }

    vector trapezium(vector times, vector func, real y0) {
        int n_step = size(times);
        vector[n_step] out;
        out[1]=y0;
        for(i in 2:n_step) {
            out[i]=out[i-1] + (times[i]-times[i-1])*(func[i-1]+func[i])/2;
        }
        return out;
    }

        real linear_primitive(real t, real t_lower, real t_upper, real usage_lower, real usage_upper) {
        real slope = (usage_upper - usage_lower)/(t_upper - t_lower);
        real s = t - t_lower;
        real out = 0.5*slope*(s^2.0) + usage_lower * s;
        return out;
    }

    real const_primitive(real t, real t_lower, real val) {
        real s = t - t_lower;
        real out = s*val;
        return out;
    }

        //integrate linear interpolation of usage function, i.e. \int_{-1}^{t_i}{ u(t) }\,dt 
    vector integrate_usage_step(vector times, array[] real ABX_usage, array[] real usage_ts, int M){
        
        //First compute the integrals of all components
        vector[M-1] component_ints;
        for(j in 2:M) {
            real t_begin = usage_ts[j-1];
            real usage_curr = ABX_usage[j-1];
            real t_end = usage_ts[j];

            component_ints[j-1] = const_primitive(t_end, t_begin, usage_curr);
        }
        
        int N = size(times);
        vector[N] out = rep_vector(0.0, N);

        for (i in 2:N) {
            real t = times[i];
            for(j in 2:M) {
                if (t > usage_ts[j-1] && t <= usage_ts[j]){
                    out[i] += const_primitive(t, usage_ts[j-1], ABX_usage[j-1]);
                    break;
                } else if (t > usage_ts[j-1]){
                    out[i] += component_ints[j-1];                    
                } else {
                    print(" ", t, ";");
                    reject("Illegal state reached");   
                }
            }
        }
        return out;
    }

    
    //integrate linear interpolation of usage function, i.e. \int_{-1}^{t_i}{ u(t) }\,dt 
    vector integrate_usage_lin(vector times, array[] real ABX_usage, array[] real usage_ts, int M){
        
        //First compute the integrals of all components
        vector[M-1] component_ints;
        for(j in 2:M) {
            real usage_prev = ABX_usage[j-1];
            real t_prev = usage_ts[j-1];

            real t_curr = usage_ts[j];
            real usage_curr = ABX_usage[j];

            component_ints[j-1] = linear_primitive(t_curr, t_prev, t_curr, usage_prev, usage_curr);
        }
        
        int N = size(times);
        vector[N] out = rep_vector(0.0, N);

        for (i in 2:N) {
            real t = times[i];
            for(j in 2:M) {
                if (t > usage_ts[j-1] && t <= usage_ts[j]){
                    out[i] += linear_primitive(t, usage_ts[j-1], usage_ts[j], ABX_usage[j-1], ABX_usage[j]);
                    break;
                } else if (t > usage_ts[j-1]){
                    out[i] += component_ints[j-1];                    
                } else {
                    reject("Illegal state reached");   
                }
            }
        }
        return out;
    }

        //integrate constant function
    vector integrate_const(vector times) {
        int N = size(times);
        vector[N] out;

        for(i in 1:N) {
            out[i] = times[i]+1.0;
        }
        return out;
    }

    vector rev_trapezium2(vector times, vector func, array[] int index) {
        int n = size(index);
        vector[n] out;
        out[1]=0; // First event is always sampling at boundary, no lineages present.

        for (i in 2:n) {
            int a = index[i-1];
            int b = index[i];
            int m = a-b;

            vector[m] dts = times[(b+1):a]-times[b:(a-1)];
            vector[m+1] fs = func[b:a];
            vector[m] summands;

            if (min(dts) < 0) reject("Invalid time vector!");

            for (j in 1:m) {
                summands[j] = 0.5*dts[j]*(fs[j]+fs[j+1]);
            }
            out[i] = sum(summands);
        }
        return out;
    }

    vector compute_rates(vector times, vector func, vector combNs, array[] int event_types, array[] int index) {
        int n = size(index);
        vector[n] out;
        out[1] = 0; // First event is always sampling at boundary, no lineages present.

        for (i in 2:n) {
            int a = index[i-1];
            int b = index[i];
            int m = a-b;

            vector[m] dts = times[(b+1):a]-times[b:(a-1)];
            vector[m+1] fs = func[b:a];
            vector[m] summands;

            if (min(dts) < 0) reject("Invalid time vector!");

            for (j in 1:m) {
                summands[j] = 0.5*dts[j]*(fs[j]+fs[j+1]);
            }
            out[i] = -combNs[i]*sum(summands);
            if (event_types[i] == 1) out[i] += (log(combNs[i]) + log(0.5*(fs[1]+fs[2])));
        }
        return out;
    }