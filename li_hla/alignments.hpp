// The class handles reading the bam file

#ifndef _LSONG_ALIGNMENT_HEADER
#define _LSONG_ALIGNMENT_HEADER

#include "samtools-0.1.19/sam.h"
#include <map>
#include <string>
#include <assert.h>
#include <iostream>

#include "defs.h"

class Alignments
{
private:
	samfile_t *fpSam ;	
	bam1_t *b ;

	char fileName[1024] ;
	bool opened ;	
	std::map<std::string, int> chrNameToId ;
	int onlyChromId ; // ignore other chromosomes

	void Open()
	{
		fpSam = samopen( fileName, "rb", 0 ) ;
		if ( !fpSam->header )
		{
			fprintf( stderr, "Can not open %s.\n", fileName ) ;
			exit( 1 ) ;
		}

		// Collect the chromosome information
		for ( int i = 0 ; i < fpSam->header->n_targets ; ++i )		
		{
			std::string s( fpSam->header->target_name[i] ) ;
			chrNameToId[s] = i ;
		}
		opened = true ;
	}
public:
	struct _pair segments[MAX_SEG_COUNT] ;		
	unsigned int segCnt ;

	Alignments() { b = NULL ; opened = false ; onlyChromId = -1 ; }
	~Alignments() {}

	void Open( char *file )
	{
		strcpy( fileName, file ) ;
		Open() ;
	}

	void Rewind()
	{
		Close() ;
		Open() ;
	}

	void Close()
	{
		samclose( fpSam ) ;
		fpSam = NULL ;
	}

	bool IsOpened()
	{
		return opened ;
	}

	int Next()
	{
		int i ;
		int start = 0, len = 0 ;
		uint32_t *rawCigar ;

		while ( 1 )
		{
			while ( 1 )
			{
				if ( b )
					bam_destroy1( b ) ;
				b = bam_init1() ;

				if ( samread( fpSam, b ) <= 0 )
					return 0 ;
				if ( b->core.flag & 0xC )
					continue ;

				if ( onlyChromId != -1 && onlyChromId != GetChromId() )
				{
					continue ;
				}

				if ( ( b->core.flag & 0x900 ) == 0 )
					break ;
			}
			// Compute the exons segments from the reads
			segCnt = 0 ;
			start = b->core.pos ; //+ 1 ;
			rawCigar = bam1_cigar( b ) ; 
			for ( i = 0 ; i < b->core.n_cigar ; ++i )
			{
				int op = rawCigar[i] & BAM_CIGAR_MASK ;
				int num = rawCigar[i] >> BAM_CIGAR_SHIFT ;

				switch ( op )
				{
					case BAM_CMATCH:
					case BAM_CDEL:
						len += num ; break ;
					case BAM_CINS:
					case BAM_CSOFT_CLIP:
					case BAM_CHARD_CLIP:
					case BAM_CPAD:
						num = 0 ; break ;
					case BAM_CREF_SKIP:
						{
							segments[ segCnt ].a = start ;
							segments[ segCnt ].b = start + len - 1 ;
							++segCnt ;
							start = start + len + num ;
							len = 0 ;
						} break ;
					default:
						len += num ; break ;
				}
			}
			segments[ segCnt ].a = start ;
			segments[ segCnt ].b = start + len - 1 ;
			++segCnt ;

			/*for ( i = 0 ; i < segCnt ; ++i )
			  printf( "(%d %d) ", segments[i].a, segments[i].b ) ;
			  printf( "\n" ) ;*/
			
			// Check whether the mates are compatible
			int mChrId = b->core.mtid ;
			int64_t mPos = b->core.mpos ;

			if ( b->core.mtid == b->core.tid )
			{
				for ( i = 0 ; i < segCnt - 1 ; ++i )
				{
					if ( mPos >= segments[i].b && mPos <= segments[i + 1].a )
						break ;
				}
				if ( i < segCnt - 1 )
					continue ;
			}
			
			break ;
		}

		return 1 ;
	}


	int GetChromId()
	{
		return b->core.tid ; 
	}

	char* GetChromName( int tid )
	{
		return fpSam->header->target_name[ tid ] ; 
	}

	int GetChromIdFromName( const char *s )
	{
		std::string ss( s ) ;
		if ( chrNameToId.find( ss ) == chrNameToId.end() )
		{
			printf( "Unknown genome name: %s\n", s ) ;
			exit( 1 ) ;
		}
		return chrNameToId[ss] ;
	}

	int GetChromLength( int tid )
	{
		return fpSam->header->target_len[ tid ] ;
	}

	void GetMatePosition( int &chrId, int64_t &pos )
	{
		chrId = b->core.mtid ;
		pos = b->core.mpos ; //+ 1 ;
	}

	int GetRepeatPosition( int &chrId, int64_t &pos )
	{
		// Look at the CC field.
		if ( !bam_aux_get( b, "CC" ) || !bam_aux_get( b, "CP" ) )
		{
			chrId = -1 ;
			pos = -1 ;
			return 0 ;
		}
		
		std::string s( bam_aux2Z( bam_aux_get(b, "CC" ) ) ) ;
		chrId = chrNameToId[ s ] ;
		pos = bam_aux2i( bam_aux_get( b, "CP" ) ) ;// Possible error for 64bit	
		return 1 ;
	}

	bool IsReverse()
	{
		if ( b->core.flag & 0x10 )	
			return true ;
		return false ;
	}

	bool IsMateReverse()
	{
		if ( b->core.flag & 0x20 )
			return true ;
		return false ;
	}

	char *GetReadId()
	{
		return bam1_qname( b ) ;
	}

	bool IsUnique()
	{
		if ( bam_aux_get( b, "NH" ) )
		{
			if ( bam_aux2i( bam_aux_get( b, "NH" ) ) > 1 )
				return false ;
		}	
		return true ;
	}

	int GetFieldI( char *f )
	{
		if ( bam_aux_get( b, f ) )
		{
			return bam_aux2i( bam_aux_get( b, f ) ) ;
		}
		return -1 ;
	}

	char *GetFieldZ( char *f )
	{
		if ( bam_aux_get( b, f ) )
		{
			return bam_aux2Z( bam_aux_get( b, f ) ) ;
		}
		return NULL ;
	}

	// -1:minus, 0: unknown, 1:plus
	int GetStrand()
	{
		if ( segCnt == 1 )
			return 0 ;
		if ( bam_aux_get( b, "XS" ) )
		{
			if ( bam_aux2A( bam_aux_get( b, "XS" ) ) == '-' )
				return -1 ;	
			else
				return 1 ;
		}
		else
			return 0 ;
	}

	void OnlyChrom( const char *chr )
	{
		onlyChromId = GetChromIdFromName( chr ) ;
	}
} ;
#endif
