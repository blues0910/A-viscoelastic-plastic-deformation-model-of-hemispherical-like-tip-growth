#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define P 1.0
#define delta 0.5
#define nu 0.5

int main(void){
	int i,zz;
	double r=1.0;
	double norm=1.0;
	double ani=1.2;
	double dt=0.004;
	double *x,*y,*alpha1,*theta,*rho,*phi,*kappas,*kappat,*sigmas,*sigmat,*GPhi,*epss,*epst,*f,*int1,*vt,*vn;
	x = (double*)calloc(10000, sizeof(double));
	y = (double*)calloc(10000, sizeof(double));
	alpha1 = (double*)calloc(10000, sizeof(double));
	theta = (double*)calloc(10000, sizeof(double));
	rho = (double*)calloc(10000, sizeof(double));
	phi = (double*)calloc(10000, sizeof(double));
	kappas = (double*)calloc(10000, sizeof(double));
	kappat = (double*)calloc(10000, sizeof(double));
	sigmas = (double*)calloc(10000, sizeof(double));
	sigmat = (double*)calloc(10000, sizeof(double));
	GPhi = (double*)calloc(10000, sizeof(double));
	epss = (double*)calloc(10000, sizeof(double));
	epst = (double*)calloc(10000, sizeof(double));
	f = (double*)calloc(10000, sizeof(double));
	int1 = (double*)calloc(10000, sizeof(double));
	vt = (double*)calloc(10000, sizeof(double));
	vn = (double*)calloc(10000, sizeof(double));

//時間発展
	int N=40;
for(i=0;i<=N;i++){
	theta[i]=M_PI*i/(2.0*N);
	x[i]=norm*r*cos(theta[i]);
	y[i]=ani*r*sin(theta[i]);
//	printf("%.20lf %.20lf\n",x[i],y[i]);
}
//x[N]=0.0;
//y[N]=ani*r;
//for(i=0;i<=N;i++){
//	printf("%.20lf %.20lf\n",x[i],y[i]);
//}
//exit(0);
for(zz=0;zz<=5000;zz++){
	//全弧長を前もって計算
	double s1=0.0;
	for(i=0;i<N;i++){
		s1=s1+sqrt((x[i+1]-x[i])*(x[i+1]-x[i])+(y[i+1]-y[i])*(y[i+1]-y[i]));
	}
//	printf("%lf\n",s1);
//	exit(0);
	double s2=0.0;
	double integrals1=0.0;
	rho[0]=0.0;
	for(i=1;i<N;i++){
		double bx1=x[i]-x[i-1];
		double by1=y[i]-y[i-1];
		double bx2=x[i+1]-x[i];
		double by2=y[i+1]-y[i];
		double bx3=x[i+1]-x[i-1];
		double by3=y[i+1]-y[i-1];
//		alpha1[i]=acos((bx1*bx2+by1*by2)/(sqrt(bx1*bx1+by1*by1)*sqrt(bx2*bx2+by2*by2)));
		alpha1[i]=fabs(fmod(atan2(by2,bx2)+2*M_PI,2*M_PI)-fmod(atan2(by1,bx1)+2*M_PI,2*M_PI));
		double l1;
		l1=sqrt(bx3*bx3+by3*by3);
		rho[i]=M_PI-alpha1[i];
		kappas[i]=2*sin(rho[i])/l1;
		phi[i]=fabs(fmod(atan2(by3,bx3)+2*M_PI,2*M_PI)-M_PI); //OK
		kappat[i]=sin(phi[i])/fabs(x[i]);
		if(i==1){
			kappas[i-1]=2*sin(2*(M_PI/2.0-phi[i]))/l1;
			kappat[i-1]=sin(M_PI/2.0)/fabs(x[0]);
		}
		if(i==1){
			sigmas[i-1]=P/(2*delta*kappat[i-1]);
			sigmat[i-1]=P*(2-kappas[i-1]/kappat[i-1])/(2*delta*kappat[i-1]);
		}
		sigmas[i]=P/(2*delta*kappat[i]);
		sigmat[i]=P*(2-kappas[i]/kappat[i])/(2*delta*kappat[i]);
		//頂点からの距離の計算
		s2=s2+sqrt(bx1*bx1+by1*by1);
		if(s1-s2<1.0)GPhi[i]=cos((s1-s2)*M_PI/2.0)*cos((s1-s2)*M_PI/2.0);
		if(s1-s2>1.0)GPhi[i]=0.0;
		double beta=2*nu*nu-2*nu+2;
		double K=sqrt(beta*sigmas[i]*sigmas[i]+beta*sigmat[i]*sigmat[i]+(beta-6*nu)*sigmas[i]*sigmat[i]);
		double sigmae=sqrt(0.5*(sigmas[i]-sigmat[i])*(sigmas[i]-sigmat[i])+0.5*(sigmat[i])*(sigmat[i])+0.5*(sigmas[i])*(sigmas[i]));
		epss[i]=GPhi[i]*0.5*sigmae*(sigmas[i]-nu*sigmat[i])/K;
		epst[i]=GPhi[i]*0.5*sigmae*(sigmat[i]-nu*sigmas[i])/K;
		double ds=sqrt(bx1*bx1+by1*by1);
		if(i==1){
			f[i]=((epss[i-1]-kappas[i-1]*epst[i-1]/kappat[i-1])/sin(M_PI/2.0)+(epss[i]-kappas[i]*epst[i]/kappat[i])/sin(phi[i]))*(ds/2.0);
			if(s1-s2<1.0)integrals1=integrals1+f[i];
		}
		if(i>1){
			f[i]=((epss[i-1]-kappas[i-1]*epst[i-1]/kappat[i-1])/sin(phi[i-1])+(epss[i]-kappas[i]*epst[i]/kappat[i])/sin(phi[i]))*(ds/2.0);
			if(s1-s2<1.0)integrals1=integrals1+f[i];
		}
		int1[i]=integrals1;
		vt[i]=sin(phi[i])*int1[i];
		vn[i]=epst[i]/kappat[i]-cos(phi[i])*int1[i];
//		printf("%d %.10lf %.10lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,x[i],y[i],alpha1[i],kappas[i],kappat[i],phi[i],sigmas[i],sigmat[i],K,s1-s2,GPhi[i],epss[i],epst[i],vt[i],vn[i],integrals1,epss[i-1],-kappas[i-1]*epst[i-1]/kappat[i-1],1.0/sin(phi[i-1]),(epss[i-1]-kappas[i-1]*epst[i-1]/kappat[i-1]));
		printf("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,x[i],y[i],alpha1[i],kappas[i],kappat[i],phi[i],sigmas[i],sigmat[i],K,s1-s2,GPhi[i],epss[i],epst[i],vt[i],vn[i],int1[i],f[i]);
//		printf("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %.10lf\n",i,x[i],y[i],kappas[i],kappat[i],sigmas[i],sigmat[i],epss[i],epst[i],vt[i],vn[i]);
	}
	if(i==N){
		//vt[N],vn[N]を求める
		x[i+1]=-x[i-1];
		y[i+1]=y[i-1];
		double bx1=x[i]-x[i-1];
		double by1=y[i]-y[i-1];
		double bx2=x[i+1]-x[i];
		double by2=y[i+1]-y[i];
		double bx3=x[i+1]-x[i-1];
		double by3=y[i+1]-y[i-1];
//		alpha1[i]=M_PI-acos((bx1*bx2+by1*by2)/(sqrt(bx1*bx1+by1*by1)*sqrt(bx2*bx2+by2*by2)));
		alpha1[i]=fabs(fmod(atan2(by2,bx2)+2*M_PI,2*M_PI)-fmod(atan2(by1,bx1)+2*M_PI,2*M_PI));
//		printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",x[i-1],y[i-1],x[i],y[i],x[i+1],y[i+1],alpha1[i],fabs(fmod(atan2(by2,bx2)+2*M_PI,2*M_PI)),fabs(fmod(atan2(by1,bx1)+2*M_PI,2*M_PI)));
		double l1;
		l1=sqrt(bx3*bx3+by3*by3);
		rho[i]=M_PI-alpha1[i];
		kappas[i]=2*sin(rho[i])/l1;
//		kappat[i]=ani/(norm*norm);//ここが原因
		kappat[i]=kappas[i];//ここが原因
//		kappat[i]=kappat[i-1]+(kappat[i-1]-kappat[i-2])*(kappat[i-1]-kappat[i-2])/(kappat[i-2]-kappat[i-3]);//ここが原因
//		kappat[i]=sin(phi[i])/fabs(x[i]);
		sigmas[i]=P/(2*delta*kappat[i]);
		sigmat[i]=P*(2-kappas[i]/kappat[i])/(2*delta*kappat[i]);
		//頂点からの距離の計算
		GPhi[i]=cos(0.0)*cos(0.0);
		double beta=2*nu*nu-2*nu+2;
		double K=sqrt(beta*sigmas[i]*sigmas[i]+beta*sigmat[i]*sigmat[i]+(beta-6*nu)*sigmas[i]*sigmat[i]);
		double sigmae=sqrt(0.5*(sigmas[i]-sigmat[i])*(sigmas[i]-sigmat[i])+0.5*(sigmat[i])*(sigmat[i])+0.5*(sigmas[i])*(sigmas[i]));
		epss[i]=GPhi[i]*0.5*sigmae*(sigmas[i]-nu*sigmat[i])/K;
		epst[i]=GPhi[i]*0.5*sigmae*(sigmat[i]-nu*sigmas[i])/K;
		double ds=sqrt(bx1*bx1+by1*by1);
//		integrals1=integrals1;
//		f[i]=2*f[i-1]-f[i-2];
//		f[i]=f[i-1]+(f[i-1]-f[i-2])*(f[i-1]-f[i-2])/(f[i-2]-f[i-3]);//ここが原因
//		if(s1-s2<1.0)integrals1=integrals1+f[i];
//		if(s1-s2<1.0)integrals1=integrals1+((epss[i-1]-kappas[i-1]*epst[i-1]/kappat[i-1])/sin(phi[i-1])+ds*((epss[i-1]-kappas[i-1]*epst[i-1]/kappat[i-1]))/sin(phi[i-1]))*(ds/2.0);
		if(s1-s2<1.0)integrals1=integrals1+((epss[i-1]-kappas[i-1]*epst[i-1]/kappat[i-1])/sin(phi[i-1])+0.0)*(ds/2.0);
		int1[i]=integrals1;
		phi[i]=0.0;
		vt[i]=0.0;
		vn[i]=epst[i]/kappat[i]-int1[i];
//		printf("%d %.10lf %.10lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,x[i],y[i],alpha1[i],kappas[i],kappat[i],phi[i],sigmas[i],sigmat[i],K,0.0,GPhi[i],epss[i],epst[i],vt[i],vn[i],integrals1,epss[i-1],-kappas[i-1]*epst[i-1]/kappat[i-1],1.0/sin(phi[i-1]),(epss[i-1]-kappas[i-1]*epst[i-1]/kappat[i-1]));
		printf("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",i,x[i],y[i],alpha1[i],kappas[i],kappat[i],phi[i],sigmas[i],sigmat[i],K,0.0,GPhi[i],epss[i],epst[i],vt[i],vn[i],int1[i],f[i],epss[i]-kappas[i]*epst[i]/kappat[i]);
	}
	//形状変化
	x[0]=x[0];
	y[0]=y[0];
	for(i=1;i<N;i++){
		x[i]=x[i]+vt[i]*cos(phi[i])*dt+vn[i]*sin(phi[i])*dt;
		y[i]=y[i]+vt[i]*sin(phi[i])*dt+vn[i]*cos(phi[i])*dt;
	}
	x[N]=0.0;
	y[N]=y[N]+vn[N]*dt;
//	for(i=0;i<=N;i++){
//		printf("%d %lf %lf %lf\n",i,x[i],y[i],gamma);
//	}
//	printf("%d %lf %lf %lf %lf\n",i,x[N],y[N],phi[N-1],rho21);
	printf("\n\n");
}
	free(x);
	free(y);
	free(alpha1);
	free(theta);
	free(rho);
	free(phi);
	free(kappas);
	free(kappat);
	free(sigmas);
	free(sigmat);
	free(GPhi);
	free(epss);
	free(epst);
	free(int1);
	free(vt);
	free(vn);
	return 0;
}
