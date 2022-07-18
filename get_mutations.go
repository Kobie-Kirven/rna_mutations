//----------------------------------------------------
// Parsing mutations in a SAM file
// 
// Author: Kobie Kirven
//----------------------------------------------------

package main

import (
	"fmt"
	"io"
	"os"
	"log"

	"github.com/biogo/hts/bam"
)

func main(){

	// Create the BAM file parser
	r, err := os.Open("/Users/lab/Downloads/sample_md.bam")
	b, err := bam.NewReader(r,2)

	if err != nil {
		log.Fatalf("could not read bam:", err)
	}
	defer b.Close()

	b.Omit(bam.AllVariableLengthData)

	var n int
	for {
		rec, err := b.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Fatalf("error reading bam: %v", err)
		}
		fmt.Println(rec)
	}

	fmt.Println(n)

}

// func sam_cigar(string cigar) {

// 	// Get the index of the mismatch in the SAM Cigar string

// }