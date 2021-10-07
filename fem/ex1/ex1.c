#include "MathTool.h"

int main ( void )
{
	Vec3d v1, v2, v3;
	VecNd  v4, v5, v6;
	Matd T1, T2, T3;

	//3次元ベクトルのテスト
	v1.x = 1;	v1.y = 2;	v1.z = 3;
	v2.x = 4;	v2.y = 5; v2.z = 6;
	printf("*************Vec3d*************\n");
	printf("v1 = ");
	printVec3( &v1 );
	printf("v2 = ");
	printVec3( &v2 );
	sumVec3andVec3( &v1, &v2, &v3 );
	printf("v1 + v2 = ");
	printVec3( &v3 );
	subVec3andVec3( &v1, &v2, &v3 );
	printf("v1 - v2 = ");
	printVec3( &v3 );
	scalingVec3( 0.5, &v1, &v3 );
	printf("0.5 * v1 = ");
	printVec3( &v3 );
	printf("v1 . v2 = ");
	printf( "%f\n", dotVec3andVec3( &v1, &v2 ) );
	crossVec3andVec3( &v1, &v2, &v3 );
	printf("v1 x v2 = ");
	printVec3( &v3 );
	printf("| v1 | = ");
	printf( "%f\n", absVec3( &v1 ) );

	//N次元ベクトルのテスト
	initVecN( &v4 );
	initVecN( &v5 );
	initVecN( &v6 );
	setVecNDim( &v4, 5 );
	setVecNDim( &v5, 5 );
	setVecNDim( &v6, 5 );
	v4.X[ 0 ] = 1;	v4.X[ 1 ] = 2;	v4.X[ 2 ] = 3;	v4.X[ 3 ] = 4;	v4.X[ 4 ] = 5;
	v5.X[ 0 ] = 6;	v5.X[ 1 ] = 7;	v5.X[ 2 ] = 8;	v5.X[ 3 ] = 9;	v5.X[ 4 ] = 10;
	printf("*************VecNd*************\n");
	printf("v4 = ");
	printVecN( &v4 );
	printf("v5 = ");
	printVecN( &v5 );
	sumVecNandVecN( &v4, &v5, &v6 );
	printf("v4 + v5= ");
	printVecN( &v6 );
	subVecNandVecN( &v4, &v5, &v6 );
	printf("v4 - v5 = ");
	printVecN( &v6 );
	scalingVecN( 0.5, &v4, &v6 );	
	printf("0.5 * v4 = ");
	printVecN( &v6 );
	printf("v4 . v5 = ");
	printf( "%f\n", dotVecNandVecN( &v4, &v5 ) );
	printf("| v4 | = ");
	printf( "%f\n", absVecN( &v4 ) );
	releaseVecN( &v4 );
	releaseVecN( &v5 );
	releaseVecN( &v6 );

	//行列のテスト
	initMat( &T1 );
	initMat( &T2 );
	initMat( &T3 );
	setMatDim( &T1, 4, 4);
	setMatDim( &T2, 4, 4);
	T1.X[ 0]=1;		T1.X[ 1]=1;		T1.X[ 2]=2;		T1.X[ 3]=2;
	T1.X[ 4]=1;		T1.X[ 5]=2;		T1.X[ 6]=5;		T1.X[ 7]=6;
	T1.X[ 8]=1;		T1.X[ 9]=2;		T1.X[10]=7;		T1.X[11]=8;
	T1.X[12]=1;		T1.X[13]=3;		T1.X[14]=8;		T1.X[15]=8;
	T2.X[ 0]=1;		T2.X[ 1]=0;		T2.X[ 2]=0;		T2.X[ 3]=0;
	T2.X[ 4]=0;		T2.X[ 5]=1;		T2.X[ 6]=0;		T2.X[ 7]=0;
	T2.X[ 8]=0;		T2.X[ 9]=0;		T2.X[10]=1;		T2.X[11]=0;
	T2.X[12]=0;		T2.X[13]=0;		T2.X[14]=0;		T2.X[15]=1;
	printf("*************Matd*************\n");
	printf("T1 = \n");
	printMat( &T1 );
	printf("T2 = \n");
	printMat( &T2 );
	sumMatandMat( &T1, &T2, &T3 );
	printf("T1 + T2 = \n");
	printMat( &T3 );
	subMatandMat( &T1, &T2, &T3 );
	printf("T1 - T2 = \n");
	printMat( &T3 );
	multiMatandMat( &T1, &T2, &T3 );
	printf("T1 * T2 = \n");
	printMat( &T3 );
	scalingMat( 0.5, &T1, &T3 );
	printf("0.5 * T1 = \n");
	printMat( &T3 );
	printf("det( T1 ) = ");
	printf( "%f\n", detMat( &T1 ) );
	trMat( &T1, &T3 );
	printf("tr( T1 ) = \n");
	printMat( &T3 );
	invMat( &T1, &T3);
	printf("T1^( -1 ) = \n");
	printMat( &T3 );
	releaseMat( &T1 );
	releaseMat( &T2 );
	releaseMat( &T3 );

	/*printf("jikken\n");
	VecNd  v7;
	initVecN( &v7 );
	setVecNDim( &v7, 4 );
	v7.X[ 0 ] = 1;	v7.X[ 1 ] = 2;	v7.X[ 2 ] = 3;	v7.X[ 3 ] = 4;
	printMat(&T1);
	printf("\n");
	printVecN(&v7);
	printf("\n");
	multiMatandVecN(&T1, &v4, &v5);
	printVecN(&v5);
	*/
	return 0;
}
