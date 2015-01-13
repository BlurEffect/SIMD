/****************************************************************
Kevin Meergans, 23040059

Multiprocessing and Parallel Systems (MAPS)
Assignment 2 - SIMD

meergans_kevin_SIMD_MAPS_2014.cpp
*****************************************************************/

#include <fstream>			//for file output
#include <iostream>			//for console output
#include <conio.h>			//for kbhit
#include "hr_time.h"		//for stopwatches
#include <stdio.h>			//for fputs
//#include <omp.h>			//for OpenMP constructs
using namespace std;

#define MAX_ROWS 2000
#define MAX_COLS 1000
#define MAX_CHARS 6

//int data[MAX_ROWS][MAX_COLS];

// align data on a 16-byte boundary
__declspec(align(16)) int data[MAX_ROWS][MAX_COLS];

CStopWatch s1,s2,s3;

void getData(void);
void sortDataOmpFor(void);
void testData(void);
void output2StringsOmpFor_itoa(void);
void output2StringsOmpFor_myASMitoa(void);
void output2stringsOmpFor_mySIMDitoaFastest(void);
void output2stringsOmpFor_mySIMDitoa(void);
void output2stringsOmpFor_mySIMDitoaExperimental(void);
void outputTimes(void);

/*
Timings done in the lab (in Release mode)

Function										  Time	     Percentual speedup (always compared to the line above)

output2StringsOmpFor_itoa()						0.183027		 0.00%
output2StringsOmpFor_myASMitoa()				0.151068		17.46%
output2stringsOmpFor_mySIMDitoa()				0.0770819		48.98%
output2stringsOmpFor_mySIMDitoaFastest() 		0.0562158		27.07%
output2stringsOmpFor_mySIMDitoaExperimental()	0.0553115		 1.61% (sometimes faster, sometimes slower than previous version)
*/


int main(void)
{

//**************************************************
	s3.startTimer();
		getData();
	s3.stopTimer();

//**************************************************
	s1.startTimer();
		sortDataOmpFor();
	s1.stopTimer();

//**************************************************
		testData();

//**************************************************
	s2.startTimer();
		//output2StringsOmpFor_itoa();			
		//output2StringsOmpFor_myASMitoa();
		//output2stringsOmpFor_mySIMDitoa();
		output2stringsOmpFor_mySIMDitoaFastest();
		//output2stringsOmpFor_mySIMDitoaExperimental();
	s2.stopTimer();

//**************************************************
		outputTimes();

	cout << "\n\nDone.";
	while (! _kbhit());  //to hold console
}


void getData()
{
	cout << "getting data...";
	srand(123); //arbitrary random number seed
	for(int i=0; i<MAX_ROWS; i++)
		for(int j=0; j<MAX_COLS; j++)
			data[i][j] = rand(); //RAND_MAX = 32767
}


//**************************************************

void sortDataOmpFor()
{
	void split_bubble(int * a, int n);

	cout << "\nsorting data...";

		//#pragma omp parallel for
		for(int i=0; i<MAX_ROWS; i++){
			//bubble sort row i
			split_bubble(data[i], MAX_COLS);
		}
}

void split_bubble(int *a, int max)
{
	//bubble sort first half of row
	for(int n=max-1; n>=max/2; n--){
        for(int j=max/2; j<n; j++){
            if(a[j] > a[j+1]){
				//swap
                int temp = a[j];
                a[j] = a[j+1];
                a[j+1] = temp;
            }
        }
	}
	//bubble sort second half of row
	for(int n=max/2 - 1; n>=0; n--){
        for(int j=0; j<n; j++){
            if(a[j] > a[j+1]){
				//swap
                int temp = a[j];
                a[j] = a[j+1];
                a[j+1] = temp;
            }
        }
	}
	//merge the two sorted halves
	int s[MAX_COLS];
	int jlow = 0;
	int jhigh = MAX_COLS/2;
	int j = 0;
	while(jlow < MAX_COLS/2 && jhigh < MAX_COLS)
	{
		if(a[jlow] <= a[jhigh])
			s[j++] = a[jlow++];
		else
			s[j++] = a[jhigh++];
	}
	while(jlow < MAX_COLS/2) s[j++] = a[jlow++];
	while(jhigh < MAX_COLS) s[j++] = a[jhigh++];

	for(j = 0; j < MAX_COLS; j++) a[j] = s[j];
}


//**************************************************
void testData()
{
	cout << "\n\ndata[0][0]                   = " << data[0][0];				//=87 for srand(123)
	cout << "\ndata[MAX_ROWS/2][MAX_COLS/2] = " << data[MAX_ROWS/2][MAX_COLS/2];	//=16440 for srand(123)
	cout << "\ndata[MAX_ROWS-1][MAX_COLS-1] = " << data[MAX_ROWS-1][MAX_COLS-1];	//=32760 for srand(123)
}

/************************************************************************************/

//--------------------------------------------------------------------------
// Original version of the output function.
// Builds two half-file strings in parallel using a standard library function.
//--------------------------------------------------------------------------
void output2StringsOmpFor_itoa()
{
	char numString[MAX_CHARS];
	string odataf1, odataf2;
	cout << "\n\noutputting data to sodata.txt...";

	//#pragma omp parallel for private (numString) schedule (static, MAX_ROWS/2)
	for(int i=0; i<MAX_ROWS; i++){
		for(int j=0; j<MAX_COLS; j++){
			_itoa_s(data[i][j],numString,10);
			if(i < MAX_ROWS/2) {
				odataf1 += numString;
				odataf1+="\t";
			}
			else {
				odataf2 += numString;
				odataf2+="\t";
			}
		}
		if(i < MAX_ROWS/2)
			odataf1+="\n";
		else
			odataf2+="\n";
	}

	FILE * sodata;
	//sodata = fopen ("sodata.txt","w");
	fopen_s(&sodata, "sodata.txt","w");
	fputs(odataf1.c_str(), sodata);
	fputs(odataf2.c_str(), sodata);
	fclose (sodata);
}

/************************************************************************************/

//--------------------------------------------------------------------------
// Original improved version of the output function.
// Builds two half-file strings in parallel using a handcrafted assembly function
// (myASMitoa).
//--------------------------------------------------------------------------
void output2StringsOmpFor_myASMitoa()
{
	void myASMitoa(int num, char * numstr);

	char numString[MAX_CHARS];
	string odataf1, odataf2;
	cout << "\n\noutputting data to sodata.txt...";

	//#pragma omp parallel for private (numString) schedule (static, MAX_ROWS/2)

	for (int i=0; i<MAX_ROWS; i++)
	{	for (int j=0; j<MAX_COLS; j++)
		{	myASMitoa (data[i][j],numString);
			if (i < MAX_ROWS/2)
			{	odataf1 += numString;
				odataf1+="\t";
			}
			else
			{	odataf2 += numString;
				odataf2+="\t";
			}
		}
		if(i < MAX_ROWS/2)
			odataf1+="\n";
		else
			odataf2+="\n";
	}

	FILE * sodata;
	//sodata = fopen ("sodata.txt","w");
	fopen_s(&sodata, "sodata.txt","w");
	fputs(odataf1.c_str(), sodata);
	fputs(odataf2.c_str(), sodata);
	fclose (sodata);
}


void myASMitoa(int num, char * numstr)
{	__asm {
			mov		ebx,numstr	// point ebx to numstr
			mov		esi, num		// store number in esi
			cmp		esi, 0			// if number is 0
			jne		nextdigit
			mov		[ebx],48		// then simply set numstr to "0"
			mov		[ebx+1],0		// add terminating null character
			jmp		enditoa			// and end

nextdigit:	mov	eax, 66666667h		// 66666667h = 2^34 / 10

			imul		esi			// edx:eax = number * (2^34 / 10)
									// therefore edx = number * (2^2 / 10)
			sar		edx,2			// edx = edx / 2^2
									// therefore edx = number / 10 (integer division)

			lea		ecx, [edx + edx*4]	// ecx = edx * 5
			add		ecx,ecx			// ecx = edx * 10
									// therefore ecx = (number div 10)*10

			mov		eax,esi				// store original number in eax
			sub		eax,ecx				// subtract ecx to leave remainder in eax
			add		eax,48				// add 48 to make eax the digit's ASCII code
			mov		[ebx],al			// store digit character in numstr
			inc		ebx

			mov		esi,edx				// move number div 10 back into esi
			cmp		esi,0				// if number div 10 = 0, we've finished
			jnz		nextdigit

			mov		[ebx],0				// so add terminating null character

			mov		edx,numstr			// the number is in reverse order of digits
nextChar:	dec		ebx					// so we need to reverse the string
			cmp		ebx,edx
			jle		enditoa
			mov		eax,[edx]
			mov		ecx,[ebx]
			mov		[ebx],al
			mov		[edx],cl
			inc		edx
			jmp     nextChar
enditoa:
	} // end of ASM block
}

/************************************************************************************/

//--------------------------------------------------------------------------
// Builds two half-file strings in parallel using the "mySIMDitoa" function
// for integer to ascii conversion. The function was left mostly unchanged to
// see the pure speedup by using SIMD. For a fully optimised version of this 
// function, see "output2stringsOmpFor_mySIMDitoaFastest" below.
//--------------------------------------------------------------------------
void output2stringsOmpFor_mySIMDitoa()
{
	void mySIMDitoa(int* num, char* numstr);

	// Will be passed to the conversion function and will store the result of the conversion
	// (4 numbers with 5 digits plus a '\t' character after each number and a '\0' character at the end -> 25 bytes)
	__declspec(align(16)) char numString[25];

	// The two half-file strings being built.
	string odataf1, odataf2;

	cout << "\n\noutputting data to sodata.txt...";

	//#pragma omp parallel for private (numString) schedule (static, MAX_ROWS/2)
	for (int i=0; i<MAX_ROWS; i++)
	{	
		// Increase j by 4 each iteration -> four numbers processed in parallel
		for (int j=0; j<MAX_COLS; j+=4)
		{	
			mySIMDitoa (&data[i][j],numString);
			if (i < MAX_ROWS/2)
			{	
				odataf1 += numString;
			}
			else
			{	
				odataf2 += numString;
			}
		}
		if(i < MAX_ROWS/2)
			odataf1+="\n";
		else
			odataf2+="\n";
	}

	FILE * sodata;
	//sodata = fopen ("sodata.txt","w");
	fopen_s(&sodata, "sodata.txt","w");
	fputs(odataf1.c_str(), sodata);
	fputs(odataf2.c_str(), sodata);
	fclose (sodata);
}

//--------------------------------------------------------------------------
// Takes an array of four integer numbers (first parameter) and converts them to 
// their ascii representations appending '\t' characters in between and adding a null
// terminator at the end. The resulting character string is written to the
// passed in character array (second parameter). The function relies on inline assembly
// code for the conversion making heavy use of SIMD instructions.
//--------------------------------------------------------------------------
void mySIMDitoa(int* num, char* numstr)
{
	// All three variables below are used for converting digits to ascii
	int intToAsciiPart1 = 16;					
	int intToAsciiPart2 = 32;
	int intToAscii      = 48;	

	int magicNumber		= 1717986919;		// Equals 10^34/10
	
	__asm 
	{
		MOV				edi, numstr				// Store a pointer to numstr in edi
		MOV				esi, num				// Store a pointer to the original numbers in esi
			
		MOV				edx, 0					// Load the ascii value for the '\0' character into edx
		MOV				ebx, 9					// Load the ascii value for the '\t' character into ebx
			
		MOVDQA			xmm0, [esi]				// Get the original numbers into the elements of xmm0
			
		VBROADCASTSS	xmm6, [intToAscii]		// Load 48 into all elements of xmm6, used for to-ascii conversion during the first iteration of the loop below
		VBROADCASTSS	xmm1, [intToAsciiPart1]	// Load 16 into all elements of xmm1, used for to-ascii conversion during further iterations of the loop
		VBROADCASTSS	xmm7, [intToAsciiPart2]	// Load 32 into all elements of xmm7, used for to-ascii conversion during further iterations of the loop
			
		VBROADCASTSS	xmm4, [magicNumber]		// Load 10^34/10 into all elements of xmm4, needed for division by 10

		MOV				ecx, 4					// Load a loop counter into ecx

		// Start of the loop

nextdigit:	
		// Get the next digit from the number and calculate the corresponding ascii value

		VPMULHUW		xmm2, xmm0, xmm4		// Multiply the current number(xmm0) with 10^34/10(xmm4) and store th result in xmm2 while maintaining the original value
		PSRAD			xmm2, 2					// Shift (number * 10^34/10)(xmm2) to the right by 2 bits -> xmm2 now holds (number DIV 10)

		VPSLLD          xmm3, xmm2, 3			// Calculate (number DIV 10)*10 by shift and add instructions and store the result in xmm3 
		PADDD			xmm3, xmm2				// (Shift left by 3 bits to multiply by eight and add original value twice to obtain the original value multiplied by 10)
		PADDD			xmm3, xmm2

		VPSUBD			xmm3, xmm0, xmm3        // Get (number % 10), that is the last digit of the current number, by subtracting (number DIV 10)*10(xmm3) from the original number(xmm0) 

		PADDD			xmm3, xmm6				// Add a certain number(xmm6) to get the ascii code associated to the digits in xmm4. If the value in xmm4 is a regular digit the addition will
												// give the corresponding ascii code for that digit. If xmm4 is zero as there are no more digits for that number, adding xmm6 will give the
												// ascii code of the empty space character
		
		// Save the digits temporarily in a register to order and accumulate them for easier writing to memory later on

		PSLLD			xmm5, 8					// Xmm5 is used to store the accumulated digits, shift its contents by 1 byte to the left to make space for the new digits
		POR				xmm5, xmm3				// Add the new digits by oring xmm5 (with its lowest byte zeroed due to the previous shifting) with xmm3 holding the digits in the lowest byte of its elements

		// Prepare the next iteration

		MOVDQA			xmm0, xmm2				// The original number for the next iteration will be (number DIV 10)(xmm2), move it into xmm0
			
		// Update xmm6 that is used to get ascii values from digits

		PXOR			xmm6, xmm6				// Set xmm6 to all zeroes
		PCMPEQD			xmm6, xmm0				// Check if the numbers for the next iteration still hold valid digits or if all digits have already been processed (only zeroes in xmm0)		
		PANDN			xmm6, xmm1				// Use negation/and to move 16(xmm1) into all elements of xmm6 where the corresponding xmm0 elements still holds valid digits

		PADDD			xmm6, xmm7				// Add 32(xmm7) to all elements of xmm6 -> elements corresponding to the xmm0 elements with valid digits will hold 48 that will convert the digits to 
												//										   their ascii representation
												//										-> elements corresponding to the xmm0 elements with zeroes in it will hold 32 that will convert the zeroes to 
												//									       the ascii representation of the empty space character

		DEC             ecx						// Decrease the loop count
		CMP				ecx, 0					// Check the loop count for 0
		JNZ				nextdigit				// If not all four iterations have been executed, jump back to the loop start
			
		// End of the loop

		// An xmm register can only hold 16 bytes, that is 16 digits. The loop above runs 4 iterations to get the first four digits of every number.
		// The code below basically executes the loop body another time in order to get the last digit (which is actually the first digit of the number) 
		// for each number and store them in another register. When writing to memory at the end the digits will be put together properly.

		VPMULHUW		xmm2, xmm0, xmm4		// Multiply the current number(xmm0) with 10^34/10(xmm4) and store the result in xmm2
		PSRAD			xmm2, 2					// Shift (number * 10^34/10)(xmm2) to the right by 2 bits -> xmm2 now holds (number DIV 10)

		VPSLLD          xmm3, xmm2, 3			// At the end of these three instructions xmm3 holds (number DIV 10)*10
		PADDD			xmm3, xmm2
		PADDD			xmm3, xmm2

		VPSUBD			xmm3, xmm0, xmm3        // Get (number % 10), that is the last digit of the current number, by subtracting (number DIV 10)*10(xmm3) from the original number(xmm0) 

		PADDD			xmm3, xmm6				// Add a certain number(xmm6) to get the ascii code associated to the digits in xmm4. If the value in xmm4 is a regular digit the addition will
												// give the corresponding ascii code for that digit. If xmm4 is zero as there are no more digits for that number, adding xmm6 will give the
												// ascii code of the empty space character

		// Write the digits back to memory

		MOV				ecx, 4					// Set the loop count for writing to memory

		// Start of the loop
nextNumber:
		// Write a number to memory

		MOVD			eax, xmm3				// Move the first digit(xmm3) of the first number to eax
		MOV				[edi], al				// Write the first digit of the number(al) to the address edi is pointing to
		MOVD			[edi+1], xmm5			// Write the next four digits of the number to edi (writes the lowest doubleword to memory)
		MOV				[edi+5], bl				// Append a tab character to the number

		// Prepare for next iteration

		PSRLDQ			xmm3, 4					// Shift xmm3 to the right by 4 bytes to move the first digit of the next number into the lowest doubleword of xmm3 
		PSRLDQ			xmm5, 4					// Shift xmm5 to the right by 4 bytes to move the four digits of the next number into the lowest doubleword of xmm5
		LEA				edi, [edi+6]			// Update the memory address to write to

		DEC             ecx						// Decrement the loop counter
		CMP				ecx, 0					// Check if the loop count is 0, that is all 4 numbers have been written to memory
		JNZ				nextNumber				// If not all four iterations have been executed, jump back to the loop start

		// End of the loop

		// Finalise

		MOV				[edi], dl				// Append a null terminator
	}
}

/************************************************************************************/

// The two character array replacing the strings to hold half-file output data.
// Have to be declared globally, as stack overflow occurs when declared locally in the function below.
// Alternative would be to declare this dynamically on the heap.
__declspec(align(16)) char odata1New[ MAX_ROWS/2 * MAX_COLS * MAX_CHARS +  MAX_ROWS/2 + 1 ];
__declspec(align(16)) char odata2New[ MAX_ROWS/2 * MAX_COLS * MAX_CHARS +  MAX_ROWS/2 + 1 ];

//--------------------------------------------------------------------------
// Builds two half-file strings in parallel using the "mySIMDitoaFastest" function
// for integer to ascii conversion. The function below has been fully optimised and 
// differs to a greater extent from the original output function. However, the general 
// structure with two strings being built has been kept though (it would have been possible to
// rely on one big character array instead of two).
//--------------------------------------------------------------------------
void output2stringsOmpFor_mySIMDitoaFastest()
{
	//__declspec(align(16)) char* odata1New = new char[ MAX_ROWS/2 * MAX_COLS * 6 +  MAX_ROWS/2 + 1 ];
	//__declspec(align(16)) char* odata2New = new char[ MAX_ROWS/2 * MAX_COLS * 6 +  MAX_ROWS/2 + 1 ];

	void mySIMDitoaFastest(int* num, char* numstr);

	cout << "\n\noutputting data to sodata.txt...";

	/*
	Restructured double for loop at the cost of some redundancy. Moved the if/else condition out of the inner for loop 
	as it is solely depending on the outer loop variable and lead to unnecessary branching within the loop. Also, the second
	if/else to append a new line to the correct string could be eliminated that way.
	*/

	/*
	Use offsets into the output char arrays to directly write converted numbers to the proper places. This way, no copying
	of the data returned by the conversion function is required. Also, using char arrays instead of strings got rid of 
	string appends, that are way more expensive.
	*/

	//#pragma omp parallel for shared (numString) schedule (static, MAX_ROWS/2)
	for ( int i = 0; i < MAX_ROWS; ++i )
	{	
		if ( i < MAX_ROWS/2 )
		{
			for ( int j = 0; j < MAX_COLS; j += 4 )
			{	
				mySIMDitoaFastest ( &data[i][j], &odata1New[i * MAX_COLS * MAX_CHARS + j * MAX_CHARS + i] );
			}
			odata1New[i * MAX_COLS * MAX_CHARS + MAX_COLS * MAX_CHARS + i] = '\n';
		}else
		{
			for ( int j = 0; j < MAX_COLS; j += 4 )
			{	
				mySIMDitoaFastest ( &data[i][j], &odata2New[( i - MAX_ROWS/2 ) * MAX_COLS * MAX_CHARS + j * MAX_CHARS + i - MAX_ROWS/2] );
			}
			odata2New[( i - MAX_ROWS/2 ) * MAX_COLS * MAX_CHARS + MAX_COLS * MAX_CHARS + ( i - MAX_ROWS/2 )] = '\n';
		}
	}

	// Append null terminators to the char arrays
	odata1New[MAX_ROWS/2 * MAX_COLS * MAX_CHARS + MAX_ROWS/2] = '\0';
	odata2New[MAX_ROWS/2 * MAX_COLS * MAX_CHARS + MAX_ROWS/2] = '\0';
	
	FILE * sodata;
	//sodata = fopen ( "sodata.txt","w" );
	fopen_s(&sodata, "sodata.txt","w");
	fputs( odata1New, sodata );
	fputs( odata2New, sodata );
	fclose ( sodata );

	//delete[] odata1New;
	//delete[] odata2New;
}


//--------------------------------------------------------------------------
// Takes an array of four integer numbers (first parameter) and converts them to 
// their ascii representations appending '\t' characters in between numbers.
// The resulting character string is written to the passed in character array (second parameter). 
// The function relies on inline assembly code for the conversion making heavy use of SIMD instructions.
// Great parts identical to "mySIMDitoa".
//--------------------------------------------------------------------------
void mySIMDitoaFastest(int* num, char* numstr)
{

	// All three variables below are used for converting digits to ascii
	int intToAsciiPart1 = 16;					
	int intToAsciiPart2 = 32;
	int intToAscii      = 48;	

	int magicNumber		= 1717986919;		// Equals 10^34/10
	
	__asm 
	{
		MOV				edi, numstr				// Store a pointer to numstr in edi
		MOV				esi, num				// Store a pointer to the original numbers in esi
			
		MOV				edx, 9					// Load the ascii value for the '\t' character into ebx
			
		MOVDQA			xmm0, [esi]				// Get the original numbers into the elements of xmm0
			
		VBROADCASTSS	xmm6, [intToAscii]		// Load 48 into all elements of xmm6, used for to-ascii conversion during the first iteration of the loop below
		VBROADCASTSS	xmm1, [intToAsciiPart1]	// Load 16 into all elements of xmm1, used for to-ascii conversion during further iterations of the loop
		VBROADCASTSS	xmm7, [intToAsciiPart2]	// Load 32 into all elements of xmm7, used for to-ascii conversion during further iterations of the loop
			
		VBROADCASTSS	xmm4, [magicNumber]		// Load 10^34/10 into all elements of xmm4, needed for division by 10

		MOV				ecx, 4					// Load a loop counter into ecx

		// Start of the loop

nextdigit:	
		// Get the next digit from the number and calculate the corresponding ascii value

		VPMULHUW		xmm2, xmm0, xmm4		// Multiply the current number(xmm0) with 10^34/10(xmm4) and store th result in xmm2 while maintaining the original value
		PSRAD			xmm2, 2					// Shift (number * 10^34/10)(xmm2) to the right by 2 bits -> xmm2 now holds (number DIV 10)

		VPSLLD          xmm3, xmm2, 3			// Calculate (number DIV 10)*10 by shift and add instructions and store the result in xmm3 
		PADDD			xmm3, xmm2				// (Shift left by 3 bits to multiply by eight and add original value twice to obtain the original value multiplied by 10)
		PADDD			xmm3, xmm2

		VPSUBD			xmm3, xmm0, xmm3        // Get (number % 10), that is the last digit of the current number, by subtracting (number DIV 10)*10(xmm3) from the original number(xmm0) 

		PADDD			xmm3, xmm6				// Add a certain number(xmm6) to get the ascii code associated to the digits in xmm4. If the value in xmm4 is a regular digit the addition will
												// give the corresponding ascii code for that digit. If xmm4 is zero as there are no more digits for that number, adding xmm6 will give the
												// ascii code of the empty space character
		
		// Save the digits temporarily in a register to order and accumulate them for easier writing to memory later on

		PSLLD			xmm5, 8					// Xmm5 is used to store the accumulated digits, shift its contents by 1 byte to the left to make space for the new digits
		POR				xmm5, xmm3				// Add the new digits by oring xmm5 (with its lowest byte zeroed due to the previous shifting) with xmm3 holding the digits in the lowest byte of its elements

		// Prepare the next iteration

		MOVDQA			xmm0, xmm2				// The original number for the next iteration will be (number DIV 10)(xmm2), move it into xmm0
			
		// Update xmm6 that is used to get ascii values from digits

		PXOR			xmm6, xmm6				// Set xmm6 to all zeroes
		PCMPEQD			xmm6, xmm0				// Check if the numbers for the next iteration still hold valid digits or if all digits have already been processed (only zeroes in xmm0)		
		PANDN			xmm6, xmm1				// Use negation/and to move 16(xmm1) into all elements of xmm6 where the corresponding xmm0 elements still holds valid digits

		PADDD			xmm6, xmm7				// Add 32(xmm7) to all elements of xmm6 -> elements corresponding to the xmm0 elements with valid digits will hold 48 that will convert the digits to 
												//										   their ascii representation
												//										-> elements corresponding to the xmm0 elements with zeroes in it will hold 32 that will convert the zeroes to 
												//									       the ascii representation of the empty space character

		DEC             ecx						// Decrease the loop count
		CMP				ecx, 0					// Check the loop count for 0
		JNZ				nextdigit				// If not all four iterations have been executed, jump back to the loop start
			
		// End of the loop

		// An xmm register can only hold 16 bytes, that is 16 digits. The loop above runs 4 iterations to get the first four digits of every number.
		// The code below basically executes the loop body another time in order to get the last digit (which is actually the first digit of the number) 
		// for each number and store them in another register. When writing to memory at the end the digits will be put together properly.

		VPMULHUW		xmm2, xmm0, xmm4		// Multiply the current number(xmm0) with 10^34/10(xmm4) and store the result in xmm2
		PSRAD			xmm2, 2					// Shift (number * 10^34/10)(xmm2) to the right by 2 bits -> xmm2 now holds (number DIV 10)

		VPSLLD          xmm3, xmm2, 3			// At the end of these three instructions xmm3 holds (number DIV 10)*10
		PADDD			xmm3, xmm2
		PADDD			xmm3, xmm2

		VPSUBD			xmm3, xmm0, xmm3        // Get (number % 10), that is the last digit of the current number, by subtracting (number DIV 10)*10(xmm3) from the original number(xmm0) 

		PADDD			xmm3, xmm6				// Add a certain number(xmm6) to get the ascii code associated to the digits in xmm4. If the value in xmm4 is a regular digit the addition will
												// give the corresponding ascii code for that digit. If xmm4 is zero as there are no more digits for that number, adding xmm6 will give the
												// ascii code of the empty space character

		// Write the digits back to memory

		MOV				ecx, 4					// Set the loop count for writing to memory

		// Start of the loop
nextNumber:

		// Write a number to memory

		MOVD			eax, xmm3				// Move the first digit(xmm3) of the first number to eax
		MOV				[edi], al				// Write the first digit of the number(al) to the address edi is pointing to
		MOVD			[edi+1], xmm5			// Write the next four digits of the number to edi (writes the lowest doubleword to memory)
		MOV				[edi+5], dl				// Append a tab character to the number

		// Prepare for next iteration

		PSRLDQ			xmm3, 4					// Shift xmm3 to the right by 4 bytes to move the first digit of the next number into the lowest doubleword of xmm3 
		PSRLDQ			xmm5, 4					// Shift xmm5 to the right by 4 bytes to move the four digits of the next number into the lowest doubleword of xmm5
		LEA				edi, [edi+6]			// Update the memory address to write to
		
		DEC             ecx						// Decrement the loop counter
		CMP				ecx, 0					// Check if the loop count is 0, that is all 4 numbers have been written to memory
		JNZ				nextNumber				// If not all four iterations have been executed, jump back to the loop start

		// End of the loop

	}	
}

/************************************************************************************/

//--------------------------------------------------------------------------
// Builds two half-file strings in parallel using the "mySIMDitoaFastest" function
// for integer to ascii conversion. The function below has been fully optimised and 
// differs to a greater extent from the original output function. However, the general 
// structure with two strings being built has been kept though (it would have been possible to
// rely on one big character array instead of two).
//--------------------------------------------------------------------------
void output2stringsOmpFor_mySIMDitoaExperimental()
{
	//__declspec(align(16)) char* odata1New = new char[ MAX_ROWS/2 * MAX_COLS * 6 +  MAX_ROWS/2 + 1 ];
	//__declspec(align(16)) char* odata2New = new char[ MAX_ROWS/2 * MAX_COLS * 6 +  MAX_ROWS/2 + 1 ];

	void mySIMDitoaExperimental(int* num, char* numstr);

	cout << "\n\noutputting data to sodata.txt...";

	/*
	Restructured double for loop at the cost of some redundancy. Moved the if/else condition out of the inner for loop 
	as it is solely depending on the outer loop variable and lead to unnecessary branching within the loop. Also, the second
	if/else to append a new line to the correct string could be eliminated that way.
	*/

	/*
	Use offsets into the output char arrays to directly write converted numbers to the proper places. This way, no copying
	of the data returned by the conversion function is required. Also, using char arrays instead of strings got rid of 
	string appends, that are way more expensive.
	*/

	//#pragma omp parallel for private (numString) schedule (static, MAX_ROWS/2)
	for ( int i = 0; i < MAX_ROWS; ++i )
	{	
		if ( i < MAX_ROWS/2 )
		{
			for ( int j = 0; j < MAX_COLS; j += 8 )
			{	
				mySIMDitoaExperimental ( &data[i][j], &odata1New[i * MAX_COLS * MAX_CHARS + j * MAX_CHARS + i] );
			}
			odata1New[i*MAX_COLS*MAX_CHARS + MAX_COLS*MAX_CHARS + i] = '\n';
		}else
		{
			for ( int j = 0; j < MAX_COLS; j += 8 )
			{	
				mySIMDitoaExperimental ( &data[i][j], &odata2New[( i - MAX_ROWS/2 ) * MAX_COLS * MAX_CHARS + j * MAX_CHARS + i - MAX_ROWS/2] );
			}
			odata2New[( i - MAX_ROWS/2 ) * MAX_COLS * MAX_CHARS + MAX_COLS * MAX_CHARS + ( i - MAX_ROWS/2 )] = '\n';
		}
	}

	// Append null terminators to the char arrays
	odata1New[MAX_ROWS/2 * MAX_COLS * MAX_CHARS + MAX_ROWS/2] = '\0';
	odata2New[MAX_ROWS/2 * MAX_COLS * MAX_CHARS + MAX_ROWS/2] = '\0';
	
	FILE * sodata;
	//sodata = fopen ( "sodata.txt","w" );
	fopen_s(&sodata, "sodata.txt","w");
	fputs( odata1New, sodata );
	fputs( odata2New, sodata );
	fclose ( sodata );

	//delete[] odata1New;
	//delete[] odata2New;
}



//--------------------------------------------------------------------------
// Given the limited range for numbers in the data array, it is possible to process
// 8 numbers in parallel (as no number is bigger than a word). This version converts
// 8 numbers to their ascii representation using "mySIMDitoaExperimental" and also includes 
// some minor changes to the previous version. It doesn't seem to be faster, though. 
// A possible reason might be the additional memory loads required and shuffle operations needed
// to convert the numbers.
//--------------------------------------------------------------------------
void mySIMDitoaExperimental(int* num, char* numstr)
{
	// Used for loading 8 integers into a single xmm register
	__declspec(align(16)) byte shuffleMask1[16] = {0,1,4,5,8,9,12,13,255,255,255,255,255,255,255,255};
	__declspec(align(16)) byte shuffleMask2[16] = {255,255,255,255,255,255,255,255,0,1,4,5,8,9,12,13};

	// used for temporarily storing digits in registers
	__declspec(align(16)) byte shuffleMask3[16] = {0,255,255,255,2,255,255,255,4,255,255,255,6,255,255,255};
	__declspec(align(16)) byte shuffleMask4[16] = {8,255,255,255,10,255,255,255,12,255,255,255,14,255,255,255};

	__declspec(align(16)) short intToAsciiPart1Array[8] = {16, 16, 16, 16, 16, 16, 16, 16};
	__declspec(align(16)) short intToAsciiPart2Array[8] = {32, 32, 32, 32, 32, 32, 32, 32};
	__declspec(align(16)) short intToAsciiArray[8]      = {48, 48, 48, 48, 48, 48, 48, 48};
	__declspec(align(16)) short magicNumberArray[8]     = {26215,26215,26215,26215,26215,26215,26215,26215};

	__asm 
	{
		MOV				edi, numstr				// Store a pointer to numstr in edi
		MOV				esi, num				// Store a pointer to the original numbers in esi
			
		MOV				edx, 9					// Load the ascii value for the '\t' character into edx
			
		// Assemble the 8 numbers in a single xmm register

		MOVDQA			xmm0, [esi]				// Get the first original numbers into the elements of xmm0
		PSHUFB			xmm0, [shuffleMask1]	// Move all relevant bytes of xmm0 into the upper quadword of xmm0
		MOVDQA			xmm1, [esi + 16]		// Get the next four original numbers into xmm1
		PSHUFB			xmm1, [shuffleMask2]	// Move all relevant bytes of xmm1 into the lower quadword of xmm1
		POR				xmm0, xmm1				// Combine the 8 original numbers in xmm0 (1 number in each word element)

		MOVDQA			xmm6, [intToAsciiArray]	// Load 48 into all elements of xmm6, used for to-ascii conversion during the first iteration of the loop below

		MOV				ecx, 4					// Load a loop counter into ecx

		// Start of the loop

nextdigit:	

		MOVDQA		    xmm1, [magicNumberArray]	

		// Get the next digit from the number and calculate the corresponding ascii value

		VPMULHUW		xmm2, xmm0, xmm1		// Multiply the current number(xmm0) with 10^34/10(xmm4) and store th result in xmm2 while maintaining the original value
		PSRLW			xmm2, 2					// Shift (number * 10^34/10)(xmm2) to the right by 2 bits -> xmm2 now holds (number DIV 10)

		VPSLLW          xmm3, xmm2, 3			// Calculate (number DIV 10)*10 by shift and add instructions and store the result in xmm3 
		PADDW			xmm3, xmm2				// (Shift left by 3 bits to multiply by eight and add original value twice to obtain the original value multiplied by 10)
		PADDW			xmm3, xmm2

		VPSUBW			xmm3, xmm0, xmm3        // Get (number % 10), that is the last digit of the current number, by subtracting (number DIV 10)*10(xmm3) from the original number(xmm0) 

		PADDW			xmm3, xmm6				// Add a certain number(xmm6) to get the ascii code associated to the digits in xmm4. If the value in xmm4 is a regular digit the addition will
												// give the corresponding ascii code for that digit. If xmm4 is zero as there are no more digits for that number, adding xmm6 will give the
												// ascii code of the empty space character
		
		// Save the digits temporarily in a register to order and accumulate them for easier writing to memory later on

		// use xmm7 for second set of numbers

		VPSHUFB			xmm4, xmm3, [shuffleMask3]
		PSLLD			xmm5, 8					// Xmm5 is used to store the accumulated digits, shift its contents by 1 byte to the left to make space for the new digits
		POR				xmm5, xmm4				// Add the new digits by oring xmm5 (with its lowest byte zeroed due to the previous shifting) with xmm3 holding the digits in the lowest byte of its elements

		VPSHUFB			xmm4, xmm3, [shuffleMask4]
		PSLLD			xmm7, 8					// Xmm5 is used to store the accumulated digits, shift its contents by 1 byte to the left to make space for the new digits
		POR				xmm7, xmm4				// Add the new digits by oring xmm5 (with its lowest byte zeroed due to the previous shifting) with xmm3 holding the digits in the lowest byte of its elements


		// Prepare the next iteration

		MOVDQA			xmm0, xmm2				// The original number for the next iteration will be (number DIV 10)(xmm2), move it into xmm0
			
		// Update xmm6 that is used to get ascii values from digits

		PXOR			xmm6, xmm6				// Set xmm6 to all zeroes
		PCMPEQW			xmm6, xmm0				// Check if the numbers for the next iteration still hold valid digits or if all digits have already been processed (only zeroes in xmm0)		
		
		MOVDQA		    xmm1, [intToAsciiPart1Array]
		PANDN			xmm6, xmm1				// Use negation/and to move 16(xmm1) into all elements of xmm6 where the corresponding xmm0 elements still holds valid digits

		MOVDQA	        xmm1, [intToAsciiPart2Array]
		PADDW			xmm6, xmm1				// Add 32(xmm1) to all elements of xmm6 -> elements corresponding to the xmm0 elements with valid digits will hold 48 that will convert the digits to 
												//										   their ascii representation
												//										-> elements corresponding to the xmm0 elements with zeroes in it will hold 32 that will convert the zeroes to 
												//									       the ascii representation of the empty space character

		DEC             ecx						// Decrease the loop count
		CMP				ecx, 0					// Check the loop count for 0
		JNZ				nextdigit				// If not all four iterations have been executed, jump back to the loop start
			
		// End of the loop


		// An xmm register can only hold 16 bytes, that is 16 digits. The loop above runs 4 iterations to get the first four digits of every number.
		// The code below basically executes the loop body another time in order to get the last digit (which is actually the first digit of the number) 
		// for each number and store them in another register. When writing to memory at the end the digits will be put together properly.

		MOVDQA			xmm1, [magicNumberArray]	

		VPMULHUW		xmm2, xmm0, xmm1		// Multiply the current number(xmm0) with 10^34/10(xmm4) and store the result in xmm2
		PSRLW			xmm2, 2					// Shift (number * 10^34/10)(xmm2) to the right by 2 bits -> xmm2 now holds (number DIV 10)

		VPSLLW          xmm3, xmm2, 3			// At the end of these three instructions xmm3 holds (number DIV 10)*10
		PADDW			xmm3, xmm2
		PADDW			xmm3, xmm2

		VPSUBW			xmm3, xmm0, xmm3        // Get (number % 10), that is the last digit of the current number, by subtracting (number DIV 10)*10(xmm3) from the original number(xmm0) 

		PADDW			xmm3, xmm6				// Add a certain number(xmm6) to get the ascii code associated to the digits in xmm4. If the value in xmm4 is a regular digit the addition will
												// give the corresponding ascii code for that digit. If xmm4 is zero as there are no more digits for that number, adding xmm6 will give the
												// ascii code of the empty space character


		// Write the digits back to memory

		// Alternate memory output using extract instructions instead of moves.
		// Unrolled output to memory

		VPEXTRB			[edi], xmm3, 14			// Extract the first digit(xmm3) of the first number and write it to memory
		VPEXTRD			[edi+1], xmm5, 0		// Write the next four digits of the number to memory
		MOV				[edi+5], dl				// Append a tab character to the number

		VPEXTRB			[edi+6], xmm3, 12		// Extract the first digit(xmm3) of the second number and write it to memory
		VPEXTRD			[edi+7], xmm5, 1		// Write the next four digits of the number to memory
		MOV				[edi+11], dl			// Append a tab character to the number

		VPEXTRB			[edi+12], xmm3, 10		// Extract the first digit(xmm3) of the third number and write it to memory
		VPEXTRD			[edi+13], xmm5, 2		// Write the next four digits of the number to memory
		MOV				[edi+17], dl			// Append a tab character to the number

		VPEXTRB			[edi+18], xmm3, 8		// Extract the first digit(xmm3) of the fourth number and write it to memory
		VPEXTRD			[edi+19], xmm5, 3		// Write the next four digits of the number to memory
		MOV				[edi+23], dl			// Append a tab character to the number

		VPEXTRB			[edi+24], xmm3, 6		// Extract the first digit(xmm3) of the fifth number and write it to memory
		VPEXTRD			[edi+25], xmm7, 0		// Write the next four digits of the number to memory
		MOV				[edi+29], dl			// Append a tab character to the number

		VPEXTRB			[edi+30], xmm3, 4		// Extract the first digit(xmm3) of the sixth number and write it to memory
		VPEXTRD			[edi+31], xmm7, 1		// Write the next four digits of the number to memory
		MOV				[edi+35], dl			// Append a tab character to the number

		VPEXTRB			[edi+36], xmm3, 2		// Extract the first digit(xmm3) of the seventh number and write it to memory
		VPEXTRD			[edi+37], xmm7, 2		// Write the next four digits of the number to memory
		MOV				[edi+41], dl			// Append a tab character to the number

		VPEXTRB			[edi+42], xmm3, 0		// Extract the first digit(xmm3) of the eigth number and write it to memory
		VPEXTRD			[edi+43], xmm7, 3		// Write the next four digits of the number to memory
		MOV				[edi+47], dl			// Append a tab character to the number

	}
}

/************************************************************************************/

// use words instead of double words???
// try to allocate big char array dynamically -> impact on timings
// get rid of last if/else in fastes Output -> 2 omp for loops
// loop to do whole row in simd (fastest version)

void outputTimes()
{
	cout << "\n\nTime for getting: " << s3.getElapsedTime();
	cout << "\n\nTime for sorting: " << s1.getElapsedTime();
	cout << "\nTime for output : " << s2.getElapsedTime();
	cout << "\nCombined time   : " << s1.getElapsedTime() + s2.getElapsedTime() + s3.getElapsedTime();
}

