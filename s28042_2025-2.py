"""
Podstawowy skrypt do łączenia się z NCBI i pobierania rekordów sekwencji genetycznych dla danego identyfikatora taksonomicznego.
"""

from io import StringIO
from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        """Initialize with NCBI credentials."""
        self.email = email
        self.api_key = api_key
        
        # Ustawienia Entrez
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'
        
    def search_taxid(self, taxid):
        """Search for all records associated with a taxonomic ID."""
        print(f"Searching for records with taxID: {taxid}")
        try:
            # Najpierw pobierz informacje taksonomiczne
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            organism_name = records[0]["ScientificName"]
            print(f"Organism: {organism_name} (TaxID: {taxid})")
            
            # Szukaj rekordów
            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y")
            search_results = Entrez.read(handle)
            count = int(search_results["Count"])
            
            if count == 0:
                print(f"No records found for {organism_name}")
                return None
                
            print(f"Found {count} records")
            
            # Zapisz wyniki wyszukiwania do późniejszego wykorzystania
            self.webenv = search_results["WebEnv"]
            self.query_key = search_results["QueryKey"]
            self.count = count
            
            return count
            
        except Exception as e:
            print(f"Error searching TaxID {taxid}: {e}")
            return None
            
    def fetch_records(self, start=0, max_records=10):
        """Fetch a batch of records using the stored search results."""
        if not hasattr(self, 'webenv') or not hasattr(self, 'query_key'):
            print("No search results to fetch. Run search_taxid() first.")
            return []
            
        try:
            # Limit, aby zapobiec przeciążeniu serwera
            batch_size = min(max_records, 500)
            
            handle = Entrez.efetch(
                db="nucleotide",
                rettype="gb",
                retmode="text",
                retstart=start,
                retmax=batch_size,
                webenv=self.webenv,
                query_key=self.query_key
            )
            
            # Surowy rekord GenBank
            records_text = handle.read()
            
            return records_text
            
        except Exception as e:
            print(f"Error fetching records: {e}")
            return ""
        
def generate_csv(records, taxid):
    accession_numbers = []
    seq_lengths = []
    descriptions = []

    for record in records:
        accession_numbers.append(record.id)
        seq_lengths.append(len(record.seq))
        descriptions.append(descriptions)

    df = pd.DataFrame({
        "accession number": accession_numbers,
        "sequence length": seq_lengths,
        "description": descriptions
    })
    path = f"taxid_{taxid}_attributes.csv"
    df.to_csv(path)
    print(f"Saved csv report to {path}.")


def generate_chart(records, taxid):
    records = sorted(list(records), key=lambda r: len(r.seq))
    accession_numbers = []
    seq_lengths = []

    for record in records:
        accession_numbers.append(record.id)
        seq_lengths.append(len(record.seq))

    fig, ax = plt.subplots()
    ax.plot(accession_numbers, seq_lengths, marker="o")
    ax.set_xlabel("Accession number")
    ax.tick_params("x", labelsize=8)
    ax.set_ylabel("Seq length")
    path = f"taxid_{taxid}_chart.png"
    plt.tight_layout()
    plt.savefig(path)
    print(f"Saved chart to {path}.")
    

def main():
    # Uzyskaj dane uwierzytelniające
    email = "s28042@pjwstk.edu.pl" # input("Enter your email address for NCBI: ")
    api_key = "9c54136959ab49f748460e18577568894f07" # input("Enter your NCBI API key: ")
    
    # Utwórz obiekt retriever
    retriever = NCBIRetriever(email, api_key)
    
    # Uzyskaj taxid od użytkownika
    taxid = "9606" # input("Enter taxonomic ID (taxid) of the organism: ")

    # [DODATEK] Pobieranie minimalnej i maksymalnej długości sekwencji
    min_len_input = input("Enter minimum sequence length (leave blank for no minimum): ").strip()
    max_len_input = input("Enter maximum sequence length (leave blank for no maximum): ").strip()
    try:
        min_len = int(min_len_input)
    except ValueError:
        min_len = None

    try:
        max_len = int(max_len_input)
    except ValueError:
        max_len = None
    
    # Szukaj rekordów
    count = retriever.search_taxid(taxid)
    
    if not count:
        print("No records found. Exiting.")
        return
        
    # Pobierz kilka pierwszych rekordów jako próbkę
    print("\nFetching sample records...")
    sample_records = retriever.fetch_records(start=0, max_records=5)

    # [DODATEK] Filtrowanie rekordów
    sio = StringIO(sample_records) # type: ignore
    records = SeqIO.parse(sio, "genbank") 
    filtered_records = []
    for record in records:
        if min_len is not None and len(record.seq) < min_len:
            continue

        if max_len is not None and len(record.seq) > max_len:
            continue

        filtered_records.append(record) 


    # [DODATEK] Generowanie csv
    generate_csv(filtered_records, taxid)

    # [DODATEK] Generowanie wykresu
    generate_chart(filtered_records, taxid)
    
    # Zapisz do pliku
    output_file = f"taxid_{taxid}_sample.gb"
    with open(output_file, "w") as f:
        f.write(sample_records)
        
    print(f"Saved sample records to {output_file}")
    print("\nNote: This is just a basic retriever. You need to extend its functionality!")

if __name__ == "__main__":
    main()