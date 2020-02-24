def time_to_jd(path,filename):
    files=sorted(glob(os.path.join(dir,filename)))
    nof=np.zeros(len(files))
    for i in range(0,len(files)):
        data=fits.open(files[i])
        header=data[0].header
        image=data[0].data
        k=np.shape(image)
        nof[i]=k[0]

        check_header=header['ACQMODE']
        
        if (check_header=='Single Scan'):
            jd_up=image
            time=header['DATE']
            t=Time(time,format='isot',scale='utc')
            time_jd=t.jd
            header.insert(15,('JD',time_jd))
            files[i]
            mod_file_1=files[i].replace('.fits','')
            fits.writeto(mod_file_1+'_sliced_'+'.fits',jd_up,header,overwrite=True)
            #print(files[i],t.jd,t.mjd,'single scan image')
        
        
        
        elif (check_header=='Kinetics'):
            exposure=header['EXPOSURE']
            print('kinetic mode image with no. of files:',files[i])
            
            name_of_file=files[i]
            mod_file=name_of_file.replace('.fits','')
            time=header['DATE']
            #print(time)
            t=Time(time,format='isot',scale='utc')
            tim=t.jd
            temp=int(nof[i])
            mod_jd=np.zeros(temp)
            exp_time=header['EXPOSURE']
            exp_time=exp_time/86400  # for the 'day' from seconds calculation.
            mod_jd[0]=tim
            for j in range(1,temp):
                mod_jd[j]=mod_jd[j-1]+exp_time
                
            for k in range(0,len(mod_jd)):
                sliced_image=image[k]
                time_jd=mod_jd[k]
                header.insert(15,('JD',time_jd))
                fits.writeto(mod_file+'_sliced_%g'%k+'.fits',sliced_image,header,overwrite=True)
                print(mod_file+'_sliced_%g'%k+'.fits has been written')

            
            
