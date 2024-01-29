import random
import os

# Fonction pour générer une séquence d'ADN aléatoire
def generate_adn_seq(length):
    # Les bases possibles de l'ADN
    bases = ['A', 'C', 'G', 'T']
    
    # pour sélectionner aléatoirement une base parmi les quatre
    return ''.join(random.choice(bases) for _ in range(length))


# Fonction pour lire une séquence d'ADN à partir d'un fichier
def lire_adn_de_fichier(nom_fichier):
    # Ouvre le fichier spécifié en mode lecture ('r')
    with open(nom_fichier, 'r') as file:
        # Lire le contenu du fichier,strip pour supprimer les espaces et sauts de ligne
        return file.read().strip()


# Function to verify the validity of the DNA sequence
def verifier_validite_adn(adn):
    valid_bases = set('ACGT')
    return all(base in valid_bases for base in adn)

# Fonction pour calculer les fréquences des bases nucléiques dans la séquence d'ADN
def calculer_frequences_bases(adn):
    # Vérifie la validité de la séquence d'ADN
    if verifier_validite_adn(adn):
        # initialiser le dictionnaire pour stocker les fréquences de chaque base
        frequencies = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

        # parcourir tous les bases
        for base in adn:
            frequencies[base] += 1

        return frequencies
    else:
        print("\nErreur : La séquence d'ADN n'est pas valide.")
        return None


# Fonction pour transcrire la séquence d'ADN en séquence d'ARN
def transcrire_adn_en_arn(adn):
    # Vérifie la validité de la séquence d'ADN
    if verifier_validite_adn(adn):
        # Remplacer 'T' par 'U' 
        return adn.replace('T', 'U')
    else:
        print("\nErreur : La séquence d'ADN n'est pas valide.")
        return None




# Fonction pour transcrire une séquence d'ARN en une séquence de protéines
def transcrire_arn_en_proteines(arn_sequence):
    # Code génétique 
    genetic_code = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    protein_sequence = ""
    # Parcours de la séquence d'ARN en sautant de 3 en 3 pour former des codons
    for i in range(0, len(arn_sequence), 3):
        codon = arn_sequence[i:i + 3]
        amino_acid = genetic_code.get(codon, '')
        # Arret de la traduction à un codon Stop
        if amino_acid == 'Stop':
            break
        protein_sequence += amino_acid

    return protein_sequence



# Fonction pour calculer le complément inverse de la séquence d'ADN
def calculer_complement_inverse_adn(adn):
    # Vérifier la validité de la séquence d'ADN
    if verifier_validite_adn(adn):
        # Créer un dictionnaire pour les compléments
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        
        # Inverser la séquence d'ADN et obtenir le complément de chaque base
        reverse_complement = ''.join(complement[base] for base in reversed(adn))
        
        return reverse_complement
    else:
        print("\nErreur: La chaîne ADN n'est pas valide.")
        return None


# Fonction pour calculer le contenu GC de la séquence d'ADN
def calculer_taux_gc(adn):
    # Vérifier la validité de la séquence d'ADN
    if verifier_validite_adn(adn):
        # Calculer le contenu GC
        nb_gc = adn.count('G') + adn.count('C')
        total_bases = len(adn)
        contenu_gc = (nb_gc / total_bases) * 100
        return contenu_gc
    else:
        print("\nErreur: La chaîne ADN n'est pas valide.")
        return None

# Fonction pour calculer les fréquences des codons dans la séquence d'ADN
def calculer_frequences_codons(adn):
    # Vérifier la validité de la séquence d'ADN
    if verifier_validite_adn(adn):
        # Initialiser un dictionnaire pour stocker les fréquences de chaque codon
        frequences = {}

        # Itérer la séquence d'ADN par blocs de 3 (codons)
        for i in range(0, len(adn) - 2, 3):
            codon = adn[i:i + 3]
            # Incrémenter la fréquence correspondante dans le dictionnaire
            frequences[codon] = frequences.get(codon, 0) + 1

        # Renvoyer le dictionnaire contenant les fréquences de chaque codon
        return frequences
    else:
        print("\nErreur: La chaîne ADN n'est pas valide.")
        return None



# Fonction pour effectuer des mutations 
def effectuer_mutations(adn, nombre_mutations):
    # Vérifier la validité de la séquence d'ADN
    if verifier_validite_adn(adn):
        # Générer des indices aléatoires pour les positions de mutation
        positions_mutation = random.sample(range(len(adn)), nombre_mutations)

        # Créer une liste à partir de la séquence d'ADN pour permettre la mutation
        adn_liste = list(adn)

        # Effectuer des mutations aux positions aléatoires
        for position in positions_mutation:
            # Obtenir la base actuelle à la position de mutation
            base_actuelle = adn_liste[position]
            
            # Générer une base de remplacement aléatoire
            base_remplacement = random.choice([base for base in 'ACGT' if base != base_actuelle])
            
            # Effectuer la mutation de substitution
            adn_liste[position] = base_remplacement

        # Convertir la liste en une chaîne de caractères
        adn_muté = ''.join(adn_liste)
        return adn_muté
    else:
        print("\nErreur: La chaîne ADN n'est pas valide.")
        return None




# Fonction pour rechercher toutes les positions d'un motif dans la séquence d'ADN
def chercher_motif(adn, motif):
    positions = []
    index = adn.find(motif)
    while index != -1:
        positions.append(index + 1)  # Ajustement pour l'indexation à partir de 1
        index = adn.find(motif, index + 1)
    return positions




# Fonction pour générer une chaîne consensus et une matrice profil 
def generer_consensus_et_profile(collection):
    # Vérifier la validité des séquences d'ADN dans la collection
    if all(verifier_validite_adn(adn) for adn in collection):
        # Obtenir la longueur des séquences d'ADN
        longueur_sequence = len(collection[0])

        # Initialiser la matrice profil
        matrice_profil = {'A': [0] * longueur_sequence, 'C': [0] * longueur_sequence, 'G': [0] * longueur_sequence, 'T': [0] * longueur_sequence}

        # Remplir la matrice profil
        for adn in collection:
            for i, base in enumerate(adn):
                matrice_profil[base][i] += 1

        # Générer la chaîne consensus
        chaine_consensus = ''.join(max(matrice_profil, key=lambda base: matrice_profil[base][i]) for i in range(longueur_sequence))

        return chaine_consensus, matrice_profil
    else:
        print("\nErreur: Une ou plusieurs séquences d'ADN de la collection ne sont pas valides.")
        return None, None



# Main program interface
while True:
    print("\nMenu:")
    print("1. Générer une chaîne ADN aléatoire")
    print("2. Choisir une chaîne ADN à partir d'un fichier")
    print("3. Quitter")

    choice = int(input("\nFaites votre choix (1, 2 ou 3): "))
    
    if choice == 1:
        length = int(input("Entrez la longueur de la chaîne ADN aléatoire: "))
        adn = generate_adn_seq(length)
        print("\nChaîne ADN aléatoire générée:", adn)
    
    elif choice == 2:
        nom_fichier = input("Entrez le nom du fichier contenant la chaîne ADN: ")
        if os.path.exists(nom_fichier):
            adn = lire_adn_de_fichier(nom_fichier)
            print("\nChaîne ADN lue à partir du fichier:", adn)
        else:
            print("Le fichier spécifié n'existe pas.")
            continue  # Revenir au début de la boucle 
    
    elif choice == 3:
        print("Goodbye!")
        break
    
    else:
        print("Choix invalide. Veuillez choisir 1, 2 ou 3.")
        continue  # Revenir au début de la boucle 
    
    # Deuxième menu
    while True:
        print("\nMenu des Traitements:")
        print("1. Vérifier la validité de la chaîne ADN")
        print("2. Calculer les fréquences des bases nucléiques")
        print("3. Transcrire la chaîne ADN en une chaîne ARN")
        print("4. Transcrire la chaîne ARN résultante en protéines")
        print("5. Calculer le complément inverse de la chaîne ADN")
        print("6. Calculer le taux de GC de la séquence ADN")
        print("7. Calculer les fréquences de codons dans la chaîne ADN")
        print("8. Réaliser des mutations aléatoires sur la chaîne ADN")
        print("9. Chercher un motif dans la chaîne ADN")
        print("10. Générer la chaîne ADN consensus et la matrice profil")


        print("9. Retour au Menu Principal")
        choix_traitement = input("Faites votre choix (1, 2, ..., 9): ")
        
        if choix_traitement == '1':
            # Verify the validity of the DNA sequence
            if verifier_validite_adn(adn):
                print("\nLa chaîne ADN est valide.")
            else:
                print("\nErreur: La chaîne ADN n'est pas valide.")
        elif choix_traitement == '2':
            # Calculer les frequences
            frequencies = calculer_frequences_bases(adn)
            if frequencies is not None:
                print("\n la chaine adn", adn)
                print("\nFréquences des bases nucléiques dans la chaîne ADN:")
                for base, frequency in frequencies.items():
                    print(f"{base}: {frequency}")
        elif choix_traitement == '3':
            # Transcrire l'adn en arn
            rna = transcrire_adn_en_arn(adn)
            if rna is not None:
             print("\nChaine adn:", adn)   
             print("\nChaîne ARN transcrit à partir de la chaîne ADN:", rna)


        elif choix_traitement == '4':
            # Transcrire l'arn en proteines
            protein_sequence = transcrire_arn_en_proteines(rna)
            print("\nChaine adn:", adn)   
            print("\nChaîne protéine résultante:", protein_sequence)
            
        elif choix_traitement == '5':
            # Calculer le complement inverse
            reverse_complement = calculer_complement_inverse_adn(adn)
            print("\nLa chaîne ADN:", adn)
            print("\nComplément inverse de la chaîne ADN:", reverse_complement)
            
        elif choix_traitement == '6':
            gc_content = calculer_taux_gc(adn)
            if gc_content is not None:
                print("\nChaine adn:", adn)   
                print("\nTaux de GC de la séquence ADN:", gc_content, "%")
        
        elif choix_traitement == '7':
            # Calculer les frequences des codons
            codon_frequencies = calculer_frequences_codons(adn)
            if codon_frequencies is not None:
                print("\nFréquences de codons dans la chaîne ADN:")
                for codon, frequency in codon_frequencies.items():
                    print(f"{codon}: {frequency}")
           
        
        elif choix_traitement == '8':
            nombre_mutations = int(input("Entrez le nombre de mutations à réaliser: "))
            mutated_adn = effectuer_mutations(adn, nombre_mutations)
            if mutated_adn is not None:
                print("\nChaine adn:", adn)  
                print(f"\nMutations aléatoires réalisées. Chaîne ADN mutée:\n{mutated_adn}")
                adn = mutated_adn  # Update the original sequence with the mutated sequence
                
        elif choix_traitement == '9':
            # chercher un motif
            motif = input("Entrez le motif à rechercher dans la chaîne ADN: ")
            motif_locations = chercher_motif(adn, motif)
            if motif_locations:
                print("\nChaine adn:", adn)  
                print(f"\nLe motif '{motif}' a été trouvé aux positions:", motif_locations)
            else:
                print(f"\nLe motif '{motif}' n'a pas été trouvé dans la chaîne ADN.")
        
        elif choix_traitement == '10':
            # Création d'une collection de séquences d'ADN pour l'analyse
            collection = [adn]
            for _ in range(9):
                 # Ajouter de 9 séquences d'ADN aléatoires supplémentaires
                collection.append(generate_adn_seq(len(adn)))  
            consensus, profile = generer_consensus_et_profile(collection)
            if consensus and profile:
                print("\nChaîne ADN Consensus:", consensus)
                print("\nMatrice Profil:")
                
                # Affichage des fréquences de chaque base pour chaque position dans la matrice profil
                for base, counts in profile.items():
                    print(f"{base}: {' '.join(map(str, counts))}")
        
        
        elif choix_traitement == '11':
            break  # Retourner au Menu Principal
        else:
            print("Choix invalide. Veuillez choisir une option valide.")
