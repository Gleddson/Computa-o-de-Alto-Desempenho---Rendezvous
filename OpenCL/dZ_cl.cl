float dZ_f(int t, float X, float gama, float ve, float H, float I, float gama_wpow, float w) {

    float wt = w*t;
    float gamat = gama*t;

    float result_J = 0.0f;
    float result1 = H * cos(wt) + I * sin(wt);
    float result2 = 0.0f;
    int n;

    for (n = 1; n <= 20; n++) {
        result_J = ve / (n * pow(X,n) * w) / (1.0f+(n * n * gama_wpow));          
			    	
		if (n%2 == 0) {                                             
			result_J = -result_J;                                   
		}                                                           
			                                                          
		result2 += result_J * pow(2.71828f, -(n * gama * t));
    }

    return (result1 - result2);
}

float I(float zl0, float gama, float X, float ve, float gama_wpow, float w) {
	
	float result_I = zl0/w - (ve/w)*(log((X+1)/X));
	float sum = 0.0f, aux = 0.0f;

	for (int n = 1; n <= 20; n++) {
		
		aux = ((ve)/(n*n*pow(X,n)*w))/(1+(n*n*gama_wpow));
		
		if (n%2 == 0) {
			aux = -aux;
		}
		
		sum += aux;
	}

	result_I += sum;

	return result_I; 
}

float H(float z0, float gama, float ve, float vexgama, float gama_wpow, float w) {
	
	float result_H = z0;
	float sum = 0.0f, aux = 0.0f;
	    	
	for (int n = 1; n <= 20; n++) {   
	    	
		//aux = ((ve * gama)/(pow(gama,n)*(w * w))) / (1 + (n * n * gama_wpow));
		aux = (ve * gama);
	    
	    if (n%2 == 0) {
	    	aux = -aux;
	    }

	    sum += aux;
	}
	    	
	result_H += sum;
	return result_H;
}  

__kernel void Kernel_dZ(                                                                                                        
   __global float* vector,
   __global float* return_dZ)                                        
   {
		                                                                   
		int t = get_global_id(0);

		return_dZ[0] = t;

		float z10 = vector[3];
	    float z0 = vector[6];

	    float w = vector[0];

	    int aux = 0;
	    int Xaux = 0;
	    float ve = 0.0f;
	    float X = 0.0f;
	       
		float result_I = 0.0f;
		float result_H = 0.0f;
		float result_dZ = 0.0f;
		   		
		float gama;
		float vexgama;
		float gama_wpow;
	   
		if (t <= 86400) {                                                         
			
			for(ve = 0.5; ve <= 5; ve += 0.5) {
				
				float vez = ve/3.0f;
	            
	            for(aux = -14; aux<=2; aux++){

	                gama = pow(10.0f, aux);
	                vexgama = vez * gama;
       				gama_wpow = (gama / w) * (gama / w);

	                for(Xaux = 1; Xaux <= 100; Xaux++) {
	                    
	                    X = Xaux;                 

	                    result_H = H(z0, gama, vez, vexgama, gama_wpow, w);
	                    result_I = I(z10, gama, X, vez, gama_wpow, w);
	                        
	                    result_dZ = dZ_f(t, X, gama, vez, result_H, result_I, gama_wpow, w); 
	                }
	            }
	        }                                                                                                                        
		}                                                                  
	}