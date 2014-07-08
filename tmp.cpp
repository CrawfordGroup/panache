
    std::shared_ptr<TwoBodyAOInt> testeri(GetERI(primary_,primary_,primary_,primary_));
    const double * testbuf = testeri->buffer();
    std::cout.precision(20);
    for(int i = 0; i < primary_->nshell(); i++)
    for(int j = 0; j < primary_->nshell(); j++)
    for(int k = 0; k < primary_->nshell(); k++)
    for(int l = 0; l < primary_->nshell(); l++)
    {
        int n = testeri->compute_shell(i,j,k,l);

        for(int m = 0; m < n; m++)
            std::cout << i << " " << j << " " << k << " " << l << " " << testbuf[m] << "\n";
    } 
