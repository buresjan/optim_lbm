//#include "ppolynomial.h"

// for file parsing
const int MAX_CHARS_PER_LINE = 1000;
const int MAX_TOKENS_PER_LINE = 100;
int PPolynomial::parse_input_file(const char* inputfile)
{
	// open file
	ifstream fin;
	fin.open(inputfile); // open a file
	if (!fin.good()) { printf("PPolynomial::parse_input_file::file %s does not exist.\n",inputfile); return -1; } 
	///CTRL->error("Mesh","Mesh","file %s does not exist",inputfile); // exit if file not found
	
	m_order=-1;
	m_pieces=-1;
	
	int lines = 0;
//	bool not_eof=true;
	// read each line of the file and determine number of elements, sides, nodes
	while (!fin.eof())
	{
		// read an entire line into memory
		char buf[MAX_CHARS_PER_LINE];
		fin.getline(buf, MAX_CHARS_PER_LINE);
	
		// parse the line into blank-delimited tokens
		int n = 0; // a for-loop index
    
		// array to store memory addresses of the tokens in buf
		const char* token[MAX_TOKENS_PER_LINE] = {}; // initialize to 0

		// parse the line
		token[0] = strtok(buf, " "); // first token
		if (token[0]) // zero if line is blank
		{
			for (n = 1; n < MAX_TOKENS_PER_LINE; n++)
			{
				token[n] = strtok(0, " "); // subsequent tokens
				if (!token[n]) break; // no more tokens
			}
		}
				
		// now read the file as we expected
		if (lines==0)
		{
			// we read the order of the polynomial first
			if (n!=1) { printf("PPolynomial::parse_input_file::line %d has incorrect number of tokens %d.\n",lines,n); return -1; } 
			sscanf(token[0],"%d",&m_order);
//			printf("order read %d\n",m_order);
		}
		
		if (lines==1)
		{
			// we read the order of the polynomial first
			if (n!=1) { printf("PPolynomial::parse_input_file::line %d has incorrect number of tokens %d.\n",lines,n); return -1; } 
			sscanf(token[0],"%d",&m_pieces);
//			printf("pieces read %d\n",m_pieces);
			// create m_divisions
			m_divisions = new double[m_pieces+1];
			// create m_coefs
			m_coefs = new double*[m_pieces];
			for (int i=0;i<m_pieces;i++) m_coefs[i] = new double[m_order];
			m_FLAG_allocated=true;
		}
		
	
		if (m_FLAG_allocated)
		{
			if (lines>=2 && lines<2+m_pieces+1)
			{
				sscanf(token[0],"%lf",&m_divisions[lines-2]);
//				printf("divisons: %d value %e\n",lines-2,m_divisions[lines-2]);
			}
			if (lines>=2+m_pieces+1)
			{
				if (n==m_order)
				{
					for (int j=0;j<m_order;j++)
					{
						sscanf(token[j],"%lf",&m_coefs[lines-2-m_pieces-1][j]);
//						printf("coefs: %d value %e\t",lines-2-m_pieces-1,m_coefs[lines-2-m_pieces-1][j]);
					}
//					printf("\n");
				}
			}
		}

/*		if (lines>=2 && not_eof) // count entities
		{
			if (strcmp(token[0],"EOF")==0) not_eof=false; else
			if (strcmp(token[0],"Element")==0)
			{
				if (strcmp(token[1],"tria")==0) tria++; else
				if (strcmp(token[1],"quad")==0) quad++; else
				CTRL->error("Mesh","parse_input_mdf_file_2D","unknown %dD object %s",m_dimension,token[1]);
				loc_entities[2]++; 
			} else
			if (strcmp(token[0],"Side")==0)
			{
				loc_entities[1]++;
			} else
			if (strcmp(token[0],"Node")==0)
			{
				loc_entities[0]++;
			} else
			if (strcmp(token[0],"Face")==0)
			{
				loc_entities[2]++;
			}
			
		}
*/		
		lines++;
	}
	fin.close();
	return 1;
}

// indexation
// -infty .. m_divisions[0] .. 0 .. m_divisions[1] .. 1 .. ... .. m_pieces-2 .. m_divisions[m_pieces-1] .. m_pieces-1 .. m_divisions[m_pieces] .. m_pieces-1 .. +infty
int PPolynomial::get_index(double x)
{
	// step 1: find the corresponding polynomial and its index i between 0 and m_pieces
	int iR = m_pieces;
	int iL = 0;
	if (x < m_divisions[0])
	{
//		printf("PPolynomial::get_index::warning x=%e is outside the bounds given by m_divisions[0]=%e\n",x,m_divisions[0]);
		return 0;
//		return -1;
	}
	if (x > m_divisions[m_pieces])
	{
//		printf("PPolynomial::get_index::warning x=%e is outside the bounds given by m_divisions[m_pieces]=%e\n",x,m_divisions[m_pieces]);
		return m_pieces-1;
//		return -1;
	}
	if (x <= m_divisions[0]) return 0;
	if (x >= m_divisions[m_pieces]) return m_pieces-1;
	int i=(int)((iR+iL)/2.0);
// int debug=0;
	while (!(m_divisions[i]<=x && m_divisions[i+1]>x))
	{
// 		debug++;
// 		printf("debug %d iL %d iR %d i %d\n",debug,iL,iR,i);
		if (m_divisions[i] <= x) iL=i;
		if (m_divisions[i] > x) iR=i;
		i=(int)((iR+iL)/2.0);
	}
	return i;
}

double PPolynomial::get_minx()
{
	if (!allocated()) return 0;
	return m_divisions[0];
}

double PPolynomial::get_maxx()
{
	if (!allocated()) return 0;
	return m_divisions[m_pieces];
}

double PPolynomial::get(double x)
{
	if (!allocated()) return 0;
	int i=get_index(x);
	if (i<0) return 0;
	double exp_x=1.0;
	double value=0;
//	i=m_pieces-1;
	// compute the polynomial
//	for (int j=0;j<m_order;j++)
	for (int j=m_order-1;j>=0;j--)
	{
//		printf("j %d value = %e plus %e (= coef %e x exp_x %e\n",j,value,m_coefs[i][j]*exp_x, m_coefs[i][j], exp_x);
		value += m_coefs[i][j]*exp_x;
		exp_x = exp_x*(x-m_divisions[i]);
	}
	return value;
}

double PPolynomial::getPeriodic(double x)
{
	if (!allocated()) return 0;
	double period = get_maxx() - get_minx();
	if (period <=0) return 0;
	int p = (int)( (x-get_minx()) / period );
	return get( get_minx() + x - p*period );
}

PPolynomial::PPolynomial(const char *inputfilename)
{
	m_FLAG_allocated=false;
	// this will create the required fields
	if (parse_input_file(inputfilename)<0) m_FLAG_allocated=false;
}

PPolynomial::~PPolynomial()
{
	if (m_FLAG_allocated) 
	{
		delete [] m_divisions;
		for (int i=0;i<m_pieces;i++) delete [] m_coefs[i];
		delete [] m_coefs;
	}
}



