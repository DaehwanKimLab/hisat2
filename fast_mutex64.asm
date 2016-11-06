;Copyright (c) 2016 Nigel Dyer
;   This accompanies the modified fast_mutex.h code and provides the assembler
;	for 64 bit MSC builds
;	rcx is the first parameter in the calling code and is a pointer to mLock
;	rdx is the second parameter in the calling code and is a pointer to oldLock

.code

tryLockAsm PROC
        mov eax,1
        xchg eax,[rcx]
        mov [rdx],eax
		ret
tryLockAsm ENDP

unlockAsm PROC
        mov eax,0
        xchg eax,[rcx]
		ret
unlockAsm ENDP

END 
