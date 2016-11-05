
.code
 
MultiplyBy10  PROC
 
    shl           RCX, 1
    mov           RAX, RCX
    shl           RCX, 2
    add           RAX, RCX
 
    ret
 
MultiplyBy10  ENDP
 

tryLockAsm PROC

;       int *ptrLock = &mLock;
;      __asm {
;        mov eax,1
;        mov ecx,ptrLock
;        xchg eax,[ecx]
;        mov oldLock,eax
;      }

        mov eax,1
        xchg eax,[rcx]
        mov [rdx],eax
		ret

tryLockAsm ENDP

unlockAsm PROC
;      int *ptrLock = &mLock;
;      __asm {
;        mov eax,0
;        mov ecx,ptrLock
;        xchg eax,[ecx]
;      }*

        mov eax,0
        xchg eax,[rcx]
		ret
unlockAsm ENDP

END 
