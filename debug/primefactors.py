def primfaktoren(number):
    n = number
    faktoren = []
    # Teile durch 2, solange möglich
    while n % 2 == 0:
        faktoren.append(2)
        n //= 2
    # Teile durch ungerade Zahlen ab 3
    faktor = 3
    while faktor * faktor <= n:
        while n % faktor == 0:
            faktoren.append(faktor)
            n //= faktor
        faktor += 2
    # Wenn am Ende noch ein Primfaktor übrig bleibt
    if n > 1:
        faktoren.append(n)
    return faktoren

# Beispiel:
zahl = 99474
print("Primfaktoren:", primfaktoren(zahl))