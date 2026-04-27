param([Parameter(ValueFromRemainingArguments=$true)][string[]]$Args)
$targetUrl = 'https://github.com/leanprover-community/ProofWidgets4/releases/download/v0.0.95+lean-v4.29.1/ProofWidgets4.tar.gz'
$source = 'C:\Users\20564\Desktop\Graduate\论文收集\code\lean4_formalization\.lake\packages\proofwidgets\.lake\ProofWidgets4-ghllkk.tar.gz'
$out = $null
for ($i = 0; $i -lt $Args.Length; $i++) {
  if ($Args[$i] -eq '-o' -and $i + 1 -lt $Args.Length) { $out = $Args[$i + 1] }
}
if ($Args -contains $targetUrl) {
  if (-not $out) { Write-Error 'curl wrapper: missing -o'; exit 2 }
  Copy-Item -LiteralPath $source -Destination $out -Force
  exit 0
}
& 'C:\Windows\System32\curl.exe' @Args
exit $LASTEXITCODE