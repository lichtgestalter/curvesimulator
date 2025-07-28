from colorama import Fore, Style, Back

print(f"{Fore.RED}This is red text{Style.RESET_ALL}")
print(f"{Style.DIM}{Fore.RED}This is red text{Style.RESET_ALL}")
print(f"{Style.NORMAL}{Fore.RED}This is red text{Style.RESET_ALL}")
print(f"{Style.BRIGHT}{Fore.RED}This is red text{Style.RESET_ALL}")
print(f"{Fore.RED}{Back.WHITE}This is red text{Style.RESET_ALL}")
print(f"{Fore.GREEN}This is green text{Style.RESET_ALL}")
print(f"{Fore.BLUE}This is blue text{Style.RESET_ALL}")


# Foreground Colors: Fore.
# Background Colors: Back.
# Possible Colors: BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE, RESET

# Text Styles: Style.
    # DIM, NORMAL, BRIGHT, RESET_ALL