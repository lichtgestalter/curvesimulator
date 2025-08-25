from colorama import Fore, Style, Back

# Foreground Colors: Fore.
# Background Colors: Back.
# Possible Colors: BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE, RESET

# Text Styles: Style.
    # DIM, NORMAL, BRIGHT, RESET_ALL

print(f"{Style.DIM}{Fore.RED}This is dimmed red text{Style.RESET_ALL}")
print(f"{Style.NORMAL}{Fore.RED}This is normal red text{Style.RESET_ALL}")
print(f"{Style.BRIGHT}{Fore.RED}This is bright red text{Style.RESET_ALL}")
print(f"{Fore.RED}{Back.WHITE}This is red text with white background{Style.RESET_ALL}")
print()
colors = ["BLACK", "RED", "GREEN", "YELLOW", "BLUE", "MAGENTA", "CYAN", "WHITE"]
color_names = Fore.BLACK, Fore.RED, Fore.GREEN, Fore.YELLOW, Fore.BLUE, Fore.MAGENTA, Fore.CYAN, Fore.WHITE

for c, cn in zip(colors, color_names):
    print(f"{cn}This is {c}{Style.RESET_ALL}")

