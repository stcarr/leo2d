function [dos,samples] = cheb_compute(input,nJobs,p,energy_rescale)
    % input = data_angle_20;
    % nJobs = 100;
    % p = 4000;
    poly_order = p;
    energy_shift = 0;
    cheb_width = 0.0001;
    interval_start = -10;
    interval_end = 10;
    num_samples = 15000;

    %{
    alpha_p = pi/(p+2);
    sample_width = cheb_width/energy_rescale;
    sample_spacing = (interval_end - interval_start)/(num_samples);

    sample_points = zeros(num_samples,1);
    for i = 1:num_samples
        sample_points(i) = ( (interval_start + sample_spacing*i) + energy_shift)/(energy_rescale);
    end

    %Determine "dampening coefficients" (independent of sample location)
    damp_coeff = zeros(poly_order + 1,1);
    for j = 1:(poly_order + 1)
        jd = j-1;
        damp_coeff(j) = ((1-jd/(p+2))*sin(alpha_p)*cos(jd*alpha_p)+(1/(p+2))*cos(alpha_p)*sin(jd*alpha_p))/sin(alpha_p);
    end

    %}
    
    T_array = zeros(1,poly_order);
    densities = zeros(num_samples,nJobs);
    
    %% Code from Paul Cazeaux

    M = num_samples;
    X = cos(pi/M*(0.5:M));
    sample_points = X;

    %%
    
    % parfor (k = 1:nJobs)
    for (k = 1:nJobs)
        T_array = input;
        
        %% Code from Paul Cazeaux    
        
        N = poly_order-1;
        mu = T_array;
        

        g = ((N+1-(0:N)).*cos(pi*(0:N)/(N+1)) + sin(pi*(0:N)/(N+1))*cot(pi/(N+1)))/(N+1);
        
        Y = idct([g(1)*mu(1); sqrt(2)*g(2:end)'.*mu(2:end)], M);
        Y = sqrt(M)*Y'./(pi*sqrt(1-X.^2));
        
        densities(:,k) = Y;
        
        %%
        
        %{

        %Now we determine the Chebyshev coefficients for a specific sample range in this loop
        
        for i = 1:num_samples
            cheb_coeff = zeros(poly_order,1);
            a = sample_points(i) - sample_width/2;
            b = sample_points(i) + sample_width/2;

            cheb_coeff(1) = (1.0/pi)*(acos(a) - acos(b));

            for j = 2:poly_order
                jd = j-1;
                cheb_coeff(j) = (2/pi)*(sin(jd*acos(a))-sin(jd*acos(b)))/jd;
            end

            %Finally we evaluate our polynomial using the T array
            T = 0;
            for j = 1:poly_order
                T = T + T_array(j)*cheb_coeff(j)*damp_coeff(j);
            end
            densities(i,k) = densities(i,k) + T/(nJobs*cheb_width);
            
        end
        %}
        
        
        %fprintf('%d \n',k);
        
    end

    dos = densities;
    samples = energy_rescale*sample_points;
end