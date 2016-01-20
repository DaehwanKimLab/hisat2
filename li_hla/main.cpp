//usage: a.out prefix_of_allele_information alignment.bam [-b backbone_id]
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <algorithm>

#include "alignments.hpp"


struct _compatible
{
	int weight ;
	double value ;
} ;

struct _snpInfo
{
	char type ; // d, s, i
	int position ;
	char nucleotide ;
	int length ;
} ;

struct _mdComponent
{
	char type ;
	int length ;
	int num ;	
} ;

struct _result
{
	int a, b ;
	double logLikelihood ;
} ;

std::map<std::string, int> snpNameToId ; 

std::map<std::string, int> alleleNameToId ;
std::vector<std::string> alleleIdToName ;


std::vector<struct _snpInfo> snpInfo ; // When using this, we already convert the snp name into snp id
std::vector< std::vector<int> > snpLink ; // What are the allele ids associate with the snp id
std::map<int, std::vector<int> > positionToSnp ; // map of the genomic coordinate to the snp id
std::vector< std::vector<int> > alleleSnpList ; // the list of snp ids associate with this allele
std::vector<int> alleleLength ; 
std::vector< struct _pair > alignmentCoords ;

bool CompResult( struct _result a, struct _result b )
{
	return a.logLikelihood > b.logLikelihood ;
}

void Split( const char *s, char delimit, std::vector<std::string> &fields )
{
	int i ;
	fields.clear() ;
	if ( s == NULL )
		return ;

	std::string f ;
	for ( i = 0 ; s[i] ; ++i )
	{
		if ( s[i] == delimit || s[i] == '\n' )	
		{
			fields.push_back( f ) ;
			f.clear() ;
		}
		else
			f.append( 1, s[i] ) ;
	}
	fields.push_back( f ) ;
	f.clear() ;
}

int main( int argc, char *argv[] )
{
	int i, j, k ;	
	FILE *fp ;
	char buffer[10100] ;
	std::vector<std::string> fields ;
	int binSize = 50 ;
	const char *backboneName = NULL ;

	// the compatilibity between the alignment and allele
	// The likelihood of this read from this allele
	struct _compatible **compatibility ; 
	int **snpAllele ; // whether a snp showed up in the allele.

	Alignments alignments ;
	alignments.Open( argv[2] ) ;

	for ( i = 3 ; i < argc ; ++i )
	{
		if ( !strcmp( argv[i], "-b" ) )
		{
			backboneName = argv[i + 1] ;
			alignments.OnlyChrom( backboneName ) ;
			++i ;
		}
		else
		{
			fprintf( stderr, "Unknown argument %s.\n", argv[i] ) ;
			exit( 1 ) ;
		}
	}

	// Parse the files associate with the snps
	// Firstly, read in the snp list have the information
	sprintf( buffer, "%s.snp", argv[1] ) ;
	fp = fopen( buffer, "r" ) ;	
	k = 0 ;
	while ( fgets( buffer, sizeof( buffer ), fp ) )
	{
		Split( buffer, '\t', fields ) ;
		if ( backboneName && strcmp( fields[2].c_str(), backboneName ) )
			continue ;
		snpNameToId[ fields[0] ] = k ;
		struct _snpInfo info ;
		info.type = fields[1][0] ;
		if ( info.type == 'd' )
		{
			info.position = atoi( fields[3].c_str() ) ;
			info.length = atoi( fields[4].c_str() ) ;
		}
		else if ( info.type == 'i' )
		{
			info.position = atoi( fields[3].c_str() ) ;
			info.length = strlen( fields[4].c_str() ) ;
		}
		else
		{
			info.position = atoi( fields[3].c_str()	) ; // notice that the snp file is 0-based index.
			info.length = 1 ;
			info.nucleotide = fields[4][0] ;
		}
		snpInfo.push_back( info ) ;
		std::vector< int > tmpList ;
		snpLink.push_back( tmpList ) ;
		for ( int p = 0 ; p < info.length ; ++p )
		{
			if ( info.type != 'i' || p == 0 )
			{
				if ( positionToSnp.find( info.position + p ) == positionToSnp.end() )
				{
					positionToSnp[ info.position + p] = tmpList ;
				}
				positionToSnp[ info.position + p].push_back( k ) ;
			}
		}
		++k ;
	} 
	fclose( fp ) ;
	// Read in the link file. Determine the id of alleles and the association
	// of alleles and snps.
	// TODO: obtain the length of each allele and take the length into account in the statistical model
	// Add the id for the backbound
	int backboneLength = 0 ;
	sprintf( buffer, "%s_backbone.fa", argv[1] ) ;
	fp = fopen( buffer, "r" ) ;
	/*for ( i = 1 ; buffer[i] && buffer[i] != ' ' && buffer[i] != '\n' ; ++i )
		;
	buffer[i] = '\0' ;
	std::string backboneName( buffer + 1 ) ;
	alleleNameToId[ backboneName ] = 0 ;
	alleleIdToName.push_back( backboneName ) ;*/
	bool start = false ;
	while ( fgets( buffer, sizeof( buffer ), fp ) )
	{
		if ( buffer[0] == '>' )
		{
			for ( i = 1 ; buffer[i] && buffer[i] != ' ' && buffer[i] != '\n' ; ++i )
				;
			buffer[i] = '\0' ;
			if ( !strcmp( backboneName, buffer + 1 ) )
			{
				start = true ;
			}
			else if ( start )
				break ;
		}
		if ( start && buffer[0] != '>' )
		{
			int len = strlen( buffer ) ;
			if ( buffer[len - 1 ] == '\n' )
				backboneLength += len - 1 ;
			else
				backboneLength += len ;
		}
	}
	fclose( fp ) ;
	/*k = 0 ;
	if ( k == 0 )
	{
		std::vector<int> tmpList ;
		alleleSnpList.push_back( tmpList ) ;
		alleleLength.push_back( backboneLength ) ;
	}*/


	// scanning the link file
	sprintf( buffer, "%s.link", argv[1] ) ;
	fp = fopen( buffer, "r" ) ;
	k = 0 ; 
	while ( fgets( buffer, sizeof( buffer ), fp ) )
	{
		std::vector<std::string> tmpFields ;
		Split( buffer, '\t', tmpFields ) ;
		// skip the snps from other backbones
		if ( snpNameToId.find( tmpFields[0] ) == snpNameToId.end() )
			continue ;

		int snpId = snpNameToId[ tmpFields[0] ] ;
		Split( tmpFields[1].c_str(), ' ', fields ) ;
		int size = fields.size() ;
	
		for ( i = 0 ; i < size ; ++i )
		{
			if ( alleleNameToId.find( fields[i] ) == alleleNameToId.end() )
			{
				//printf( "%s %d\n", fields[i].c_str(), k ) ;
				alleleNameToId[ fields[i] ] = k ;
				alleleIdToName.push_back( fields[i] ) ;
				std::vector<int> tmpList ;
				alleleSnpList.push_back( tmpList ) ;
				alleleLength.push_back( backboneLength ) ;
				++k ;
			}

			int alleleId = alleleNameToId[ fields[i] ] ;
			//if ( snpId == 118 )
			//	printf( "%s: %s %d\n", tmpFields[0].c_str(), fields[i].c_str(), alleleId ) ;
			snpLink[ snpId ].push_back( alleleId ) ;
			alleleSnpList[ alleleId ].push_back( snpId ) ;
			if ( snpInfo[ snpId ].type == 'd' )
			{
				alleleLength[ alleleId ] -= snpInfo[ snpId ].length ;
			}
			else if ( snpInfo[ snpId ].type == 'i' )
			{
				alleleLength[ alleleId ] += snpInfo[ snpId ].length ;
			}
		}
	} 
	fclose( fp ) ;
	int numOfAllele = alleleIdToName.size() ;
	int numOfSnps = snpLink.size() ;
	snpAllele = new int* [numOfSnps] ;
	for ( i = 0 ; i < numOfSnps ; ++i )
	{
		snpAllele[i] = new int[numOfAllele] ;
		memset( snpAllele[i], 0, sizeof( int ) * numOfAllele ) ;
	}
	for ( i = 0 ; i < numOfSnps ; ++i )
	{
		int size = snpLink[i].size() ;
		for ( j = 0 ; j < size ; ++j )
		{
			snpAllele[i][ snpLink[i][j] ] = 1 ;
		}
	}

	// Compute the compatbility score for each alignment and the allele
 	// Get the number of alignment
	int numOfAlignments = 0 ;
	while ( alignments.Next() )
		++numOfAlignments ;
		
	alignments.Rewind() ;

	compatibility = new struct _compatible*[numOfAlignments] ;	
	for ( i = 0 ; i < numOfAlignments ; ++i )
	{
		compatibility[i] = new struct _compatible[ numOfAllele ] ;
		for ( j = 0 ; j < numOfAllele ; ++j )
		{
			compatibility[i][j].value = 0 ;//-log( (double)alleleLength[j] ) / log( 10.0 );
		}
	}

	i = 0 ; 
	bool *snpHit = new bool[ numOfSnps ] ;
	while ( alignments.Next() )
	{
		struct _pair coord = alignments.segments[0] ;
		alignmentCoords.push_back( coord ) ;
		memset( snpHit, 0, sizeof( bool ) * numOfSnps ) ;
		Split( alignments.GetFieldZ( "Zs" ), ',', fields ) ;
		int size = fields.size() ;
		for ( k = 0 ; k < size ; ++k )
		{
			std::vector<std::string> subfields ;
			Split( fields[k].c_str(), '|', subfields ) ;
			int snpId = snpNameToId[ subfields[2] ] ;	
			snpHit[ snpId ] = true ;
		}
		
		for ( k = coord.a ; k <= coord.b ; ++k )
		{
			int size = positionToSnp[k].size() ;
			for ( int l = 0 ; l < size ; ++l )
			{
				// if this SNP is hit. Then other allele don't have this snp 
				// will deduct its likelihood
				//TODO: the deduction can be based on the quality score of the read
				int tag = 0 ;
				int snpId = positionToSnp[k][l] ;

				if ( snpHit[ snpId ] )
				{
					tag = 0 ;
				}
				else
				{
				// if this SNP is not hit, then every allele containing this snp 
				// will deduct its likelihood
					tag = 1 ;
				}
				for ( j = 0 ; j < numOfAllele ; ++j )
				{
					if ( snpAllele[ snpId  ][j] == tag )
					{
						int v = -2 ;
						//if ( snpInfo[ snpId ].type == 'd' || snpInfo[ snpId ].type == 'i' )
						//	v = -4 * snpInfo[ snpId ].length ;	
						if ( snpInfo[ snpId ].type == 'd' && snpInfo[ snpId ].position < k 
							&& k != coord.a )
						{
							// The penality has already been subtracted.
							v = 0 ;
						}
						compatibility[i][j].value += v ;
						/*if ( i == 8 && j == 78 )
						{
							printf( "Bad snp %d: %d %d\n", tag, k, positionToSnp[k][l] ) ;
						}*/
					}
				}
			}
		}
		++i ;
	}
	//printf( "%d %d\n", numOfAlignments, numOfAllele ) ;
	// Now, let's consider every pair of alleles, and compute its log likelihood 
	double **logLikelihood ;
	logLikelihood = new double *[ numOfAllele] ;
	for ( j = 0 ; j < numOfAllele ; ++j )
	{
		logLikelihood[j] = new double[ numOfAllele ] ;
		//memset( logLikelihood[j], 0, sizeof( double ) * numOfAllele ) ;
		for ( k = 0 ; k < numOfAllele ; ++k )
			logLikelihood[j][k] = 0 ;
	}

	int prevBin = -1 ;
	double assignJBin = 0 ;
	double assignKBin = 0 ;
	for ( j = 0 ; j < numOfAllele ; ++j )
	{
		for ( k = j ; k < numOfAllele ; ++k )
		{
			double binAdjust = 0 ;
			double averageRead = ( (double)numOfAlignments ) / (double)( alleleLength[j] + alleleLength[k] ) * binSize ;
			for ( i = 0 ; i < numOfAlignments ; ++i )
			{
				double vj = compatibility[i][j].value ;
				double vk = compatibility[i][k].value ;
				double weightJ = 0, weightK = 0 ;
				if ( vj == vk )
				{
					weightJ = weightK = 0.5 ;
				}
				else if ( vj == vk + 2 )
				{
					if ( vj == 0 )
					{
						weightJ = 1 ;
					}
					else
					{
						weightJ = 0.99 ;
						weightK = 0.01 ;
					}
				}
				else if ( vk == vj + 2 )
				{
					if ( vk == 0 )
					{
						weightK = 1 ;
					}
					else
					{
						weightJ = 0.01 ;
						weightK = 0.99 ;
					}
				}
				else
				{
					if ( vk > vj )
						weightK = 1 ;
					else
						weightJ = 1 ;
				}
				
				double l = weightJ * compatibility[i][j].value + weightK * compatibility[i][k].value ;
				if ( alignmentCoords[i].a / binSize != prevBin )
				{
					if ( prevBin != -1 && 
						( assignJBin > averageRead + 4 * sqrt( averageRead ) 
							|| assignKBin > averageRead + 4 * sqrt( averageRead ) ) )
					{
						//if ( j == 8 && k == 78 )
						//	printf( "%lf: %lf %lf %d %d\n", averageRead, assignJBin, assignKBin, alleleLength[j], alleleLength[k] ) ;
						binAdjust -= 4 ;
					}
					prevBin = alignmentCoords[i].a / binSize ;
					assignJBin = 0 ;
					assignKBin = 0 ;
				}
				assignJBin += weightJ ;
				assignKBin += weightK ;
				
				/*if ( j == 8 && k == 78 && l < 0 )
				{
					printf( "Bad alignment %d (%s %s). %lf %lf: %lf\n", i,
						alleleIdToName[j].c_str(), alleleIdToName[k].c_str(),
						 compatibility[i][j].value, compatibility[i][k].value, l ) ;
				}*/
				logLikelihood[j][k] += l ;
			}
			logLikelihood[j][k] += ( -log( (double)alleleLength[j] ) / log(10.0 ) -
							log( (double)alleleLength[k] ) / log(10.0) ) ;
			logLikelihood[j][k] += binAdjust ;
		}
	}
	
	// Find the result
	double max ;
	int maxj = -1 ;
	int maxk = -1 ;
	std::vector< struct _result > results ;
	for ( j = 0 ; j < numOfAllele ; ++j )
	{
		for ( k = j ; k < numOfAllele ; ++k )
		{
			if ( maxj == -1 || logLikelihood[j][k] > max )
			{
				maxj = j ;
				maxk = k ;
				max = logLikelihood[j][k] ;
			}
			struct _result r ;
			r.a = j ;
			r.b = k ;
			r.logLikelihood = logLikelihood[j][k] ;
			results.push_back( r ) ;
		}
	}

	//printf( "%s %s %lf\n", alleleIdToName[ maxj ].c_str(), alleleIdToName[ maxk ].c_str(), max) ;
	//printf( "%lf\n", logLikelihood[124][128] ) ;	
	
	if ( results.size() == 0 )
	{
		printf( "-1 -1 -1\n" ) ;
		exit( 1 ) ;
	}
	std::sort( results.begin(), results.end(), CompResult ) ;	
	i = 0 ;
	printf( "%s %s %lf\n", alleleIdToName[ results[i].a ].c_str(), alleleIdToName[ results[i].b ].c_str(), 
			results[i].logLikelihood ) ;
	k = results.size() ;
	for ( i = 1 ; i < k ; ++i )
	{
		if ( results[i].logLikelihood != results[0].logLikelihood )
			break ;
		printf( "%s %s %lf\n", alleleIdToName[ results[i].a ].c_str(), alleleIdToName[ results[i].b ].c_str(), 
			results[i].logLikelihood ) ;
	}
	return 0 ;
}
